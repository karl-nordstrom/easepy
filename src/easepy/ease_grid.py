import logging
import typing

import numpy as np
import pyproj

logger = logging.getLogger(__name__)


class EaseGrid(object):
    """
    EASE Grid class
    The Equal Area Scaleable Earth (EASE) grid is an equal area grid for earth
    observation data. It has three different projections:
        * Northern Hemisphere
        * Southern Hemisphere
        * Global (not defined for |lat| > 84 degrees)
    The EASE grid is defined using a two-dimensional coordinate system. We denote these
    dimensions x and y in this code. These are defined in such a way that
    a regular grid in this coordinate system produces a high-quality equal area grid within
    the region of validity. The EASE grid indices (here denoted col_ix, row_iy)
    are therefore determined as follows:
        * Convert geodetic coordinate to requested EASE coordinate system
        * Determine grid indices in the x and y coordinates using a regular grid
    Due to this procedure, there is some ambiguity in what we mean by conversion functions
    such as geodetic2ease and ease2geodetic: is the user requesting a conversion to the
    EASE coordinate system (or back), or are they working with grid indices?
    The geodetic2ease conversion returns both the grid indices and the coordinates in the
    EASE coordinate system. This is done because the common use-case is to grid values into
    EASE grids, but the EASE coordinates have to be calculated as part of this, and it is therefore
    unnecessary to maintain separate functions.
    When converting back to geodetic, two functions are instead provided:
        * ease_coord2geodetic
        * ease_index2geodetic
    These convert EASE coordinates back to geodetic (always correct, regardless of the setup of the grid,
    as long as the same projection is used) or the EASE grid index back to geodetic (only correct
    if the current grid is set up with the same projection AND resolution).
    Description and further info: https://nsidc.org/ease/ease-grid-projection-gt
    Values are assumed to be in meters and degrees.
    """

    def __init__(self, resolution_m: int, projection: str) -> None:
        """
        Parameters
        ----------
        resolution_m : int
            resolution in meters
        projection : str
            projection to be used (NorthHemi, SouthHemi, or Global)
        """
        # Casting resolution to integer since needs to be a multiple of 3
        self.resolution = int(resolution_m)
        self.projection = projection
        # lat/lon projection
        #  Note that pyproj uses a lon, lat convention in argument order
        self.proj_latlon = pyproj.Proj("EPSG:4326")

        self.easeGL3km_xcols = 11568
        self.easeGL3km_yrows = 4872

        self.easeNH3km_xcols = 6000
        self.easeNH3km_yrows = 6000

        self.easeSH3km_xcols = 6000
        self.easeSH3km_yrows = 6000

        if self.resolution % 3 != 0:
            invalid_res = self.resolution
            self.resolution = -1
            msg = f"Unsupported resolution {invalid_res} meters! (only multiples of 3 meters allowed)"
            raise ValueError(msg)

        res_scale = int(self.resolution / 3000)

        if self.projection == "NorthHemi":
            self.description = "EASE Northern Hemisphere, Lambert Azimuthal projection"
            self.proj_sgrid = pyproj.Proj("EPSG:6931")
            # The validity range of the grid in terms of EASE coordinates
            # Defined here in terms of the coordinates of one of the corners of the grid
            #  These are symmetric so both ranges will be -|min/max| < 0 < |min/max|
            self._xmin, self._ymax = (-9000000.0, 9000000.0)
            self.number_cols = int(self.easeNH3km_xcols / res_scale)
            self.number_rows = int(self.easeNH3km_yrows / res_scale)

        elif self.projection == "SouthHemi":
            self.description = "EASE Southern Hemisphere, Lambert Azimuthal projection"
            self.proj_sgrid = pyproj.Proj("EPSG:6932")
            self._xmin, self._ymax = (-9000000.0, 9000000.0)
            self.number_cols = int(self.easeSH3km_xcols / res_scale)
            self.number_rows = int(self.easeSH3km_yrows / res_scale)

        elif self.projection == "Global":
            self.description = "EASE Global, Equal-Area projection"
            self.proj_sgrid = pyproj.Proj("EPSG:6933")
            self._xmin, self._ymax = (-17367530.45, 7314540.83)
            self.number_cols = int(self.easeGL3km_xcols / res_scale)
            self.number_rows = int(self.easeGL3km_yrows / res_scale)
        else:
            msg = f"Unsupported projection {self.projection}! (must be NorthHemi/SouthHemi/Global)"
            raise ValueError(msg)

        self.map_res_x = np.abs(self._xmin * 2) / self.number_cols
        self.map_res_y = np.abs(self._ymax * 2) / self.number_rows

        logger.debug(f" {self.description}")
        logger.debug(f" Resolution: {self.resolution} m")
        logger.debug(f" Number of Columns: {self.number_cols}")
        logger.debug(f" Number of Rows: {self.number_rows}")
        logger.debug(f" Res. in the x direction: {self.map_res_x} m")
        logger.debug(f" Res. in the y direction: {self.map_res_y} m")

        # These are the actual coordinate converters
        self.trans_lonlat2xy = pyproj.Transformer.from_proj(
            self.proj_latlon, self.proj_sgrid, always_xy=True
        )

        self.trans_xy2lonlat = pyproj.Transformer.from_proj(
            self.proj_sgrid, self.proj_latlon, always_xy=True
        )

    def geodetic2ease(
        self, lat: np.ndarray, lon: np.ndarray
    ) -> typing.Tuple[
        typing.Tuple[np.ndarray, np.ndarray], typing.Tuple[np.ndarray, np.ndarray]
    ]:
        """
        Function to find corresponding EASE coordinates and grid index for given lat/lon point
        Parameters
        ----------
        lat : np.ndarray
            latitude(s) of the point(s) (in degrees)
        lon : np.ndarray
            longitude(s) of the point(s) (in degrees)
        Returns
        -------
        ease_coords : tuple[xcol_id, yrow_id], tuple[xx, yy]
            EASE grid indices of the point(s), and corresponding EASE projection coordinates.
            Take same shape as input values.
        """

        if (np.abs(np.array([lat])) > 90.0).any():
            msg = "There are input lat values with absolute values above 90 degrees."
            raise ValueError(msg)

        # find EASE projection x and y coordinates of the point at lon, lat
        #  Note lon/lat convention used by pyproj is opposite to our own
        xx, yy = self.trans_lonlat2xy.transform(lon, lat)

        # check max and min values to make sure the points are within the grid
        if np.array([(xx < self._xmin)]).any() or np.array([(xx > -self._xmin)]).any():
            msg = "Some geodetic coordinates are outside of EASE grid validity range. \
                   Check documentation at https://nsidc.org/ease/ease-grid-projection-gt."
            raise ValueError(msg)
        if np.array([(yy > self._ymax)]).any() or np.array([(yy < -self._ymax)]).any():
            msg = "Some geodetic coordinates are outside of EASE grid validity range. \
                   Check documentation at https://nsidc.org/ease/ease-grid-projection-gt."
            raise ValueError(msg)

        # find EASE grid indexes of the point
        xcol_id = ((xx - self._xmin) / self.map_res_x).astype(int)
        yrow_id = ((self._ymax - yy) / self.map_res_y).astype(int)

        return (xcol_id, yrow_id), (xx, yy)

    def ease_coord2geodetic(
        self, xx: np.ndarray, yy: np.ndarray
    ) -> typing.Tuple[np.ndarray, np.ndarray]:
        """
        Function to find corresponding lat/lon point for given EASE point.
        Valid as long as the projection used is consistent.
        Parameters
        ----------
        xx : np.ndarray
            x coordinate in EASE grid coordinate system
        yy : np.ndarray
            y coordinate in EASE grid coordinate system
        Returns
        -------
        lat, lon : tuple[np.ndarray, np.ndarray]
            lat and lon values in geodetic coordinate system
        """

        #  Note lon/lat convention used by pyproj is opposite to our own
        lon, lat = self.trans_xy2lonlat.transform(xx, yy)
        return lat, lon

    def ease_index2geodetic(
        self, xcol_id: np.ndarray, yrow_id: np.ndarray
    ) -> typing.Tuple[np.ndarray, np.ndarray]:
        """
        Function to find corresponding lat/lon point for given EASE grid index pair
        Valid only if the projection and resolution used are the same.
        Returns the location of the approximate midpoint of the grid cell.
        Parameters
        ----------
        col_ix : np.ndarray
            x coordinate in EASE grid coordinate system
        row_iy : np.ndarray
            y coordinate in EASE grid coordinate system
        Returns
        -------
        lat, lon : tuple[np.ndarray, np.ndarray]
            lat and lon values in geodetic coordinate system
        """

        xx = (xcol_id + 0.5) * self.map_res_x + self._xmin
        yy = -(yrow_id + 0.5) * self.map_res_y + self._ymax

        #  Note lon/lat convention used by pyproj is opposite to our own
        lon, lat = self.trans_xy2lonlat.transform(xx, yy)
        return lat, lon
