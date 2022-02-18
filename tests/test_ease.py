from pathlib import Path
from unittest import TestCase

import numpy as np
import pytest

from easepy import EaseGrid

path = Path(__file__)


class TestEasepy(TestCase):
    def setUp(self):
        pass

    @pytest.mark.unit
    def test_ease_north_hemi_geodetic2ease_array(self):
        ease = EaseGrid(12000, "NorthHemi")
        lats = np.array([75, 85, 89.99])
        lons = np.array([-175, 7, 155])
        (x_ind, y_ind), (x, y) = ease.geodetic2ease(lats, lons)
        self.assertTrue(
            np.isclose(
                x, np.array([-145571.88051324, 68037.02296194, 472.03909014])
            ).all()
        )
        self.assertTrue(
            np.isclose(
                y, np.array([1.66389421e06, -5.54117085e05, 1.01229124e03])
            ).all()
        )
        self.assertTrue((x_ind == np.array([737, 755, 750])).all())
        self.assertTrue((y_ind == np.array([611, 796, 749])).all())

    @pytest.mark.unit
    def test_ease_north_hemi_geodetic2ease_array_outside_grid(self):
        ease = EaseGrid(12000, "NorthHemi")
        lats = np.array([-80, 85, 89.99])
        lons = np.array([-175, 7, 155])
        with pytest.raises(Exception) as e_info:
            (x_ind, y_ind), (x, y) = ease.geodetic2ease(lats, lons)
        self.assertEqual(e_info.type, ValueError)

    @pytest.mark.unit
    def test_ease_north_hemi_geodetic2ease_array_invalid_coordinate(self):
        ease = EaseGrid(12000, "NorthHemi")
        lats = np.array([95, 85, 89.99])
        lons = np.array([-175, 7, 155])
        with pytest.raises(Exception) as e_info:
            (x_ind, y_ind), (x, y) = ease.geodetic2ease(lats, lons)
        self.assertEqual(e_info.type, ValueError)

    @pytest.mark.unit
    def test_ease_north_hemi_geodetic2ease_outside_grid(self):
        with pytest.raises(Exception) as e_info:
            ease = EaseGrid(12000, "NorthHemi")
            lats = -10
            lons = 180
            (x_ind, y_ind), (x, y) = ease.geodetic2ease(lats, lons)
        self.assertEqual(e_info.type, ValueError)

    @pytest.mark.unit
    def test_ease_south_hemi_geodetic2ease_outside_grid(self):
        with pytest.raises(Exception) as e_info:
            ease = EaseGrid(12000, "SouthHemi")
            lats = 10
            lons = 180
            (x_ind, y_ind), (x, y) = ease.geodetic2ease(lats, lons)
        self.assertEqual(e_info.type, ValueError)

    @pytest.mark.unit
    def test_ease_global_geodetic2ease_outside_grid(self):
        with pytest.raises(Exception) as e_info:
            ease = EaseGrid(12000, "SouthHemi")
            lats = 89
            lons = 180
            (x_ind, y_ind), (x, y) = ease.geodetic2ease(lats, lons)
        self.assertEqual(e_info.type, ValueError)

    @pytest.mark.unit
    def test_ease_raises_resolution_value_error(self):
        with pytest.raises(Exception) as e_info:
            ease = EaseGrid(11000, "SouthHemi")
            ease.geodetic2ease(50, 100)
        self.assertEqual(e_info.type, ValueError)

    @pytest.mark.unit
    def test_ease_raises_projection_value_error(self):
        with pytest.raises(Exception) as e_info:
            ease = EaseGrid(12000, "SH")
            ease.geodetic2ease(50, 100)
        self.assertEqual(e_info.type, ValueError)

    @pytest.mark.unit
    def test_ease_raises_latitude_value_error(self):
        with pytest.raises(Exception) as e_info:
            ease = EaseGrid(12000, "SouthHemi")
            ease.geodetic2ease(91, 100)
        self.assertEqual(e_info.type, ValueError)

    @pytest.mark.unit
    def test_ease_north_hemi_ease_coord2geodetic(self):
        ease = EaseGrid(12000, "NorthHemi")
        yy = -554117.0849300659
        xx = 68037.02296194
        lat, lon = ease.ease_coord2geodetic(xx, yy)
        self.assertTrue(np.isclose(np.array([lat]), np.array([85])).all())
        self.assertTrue(np.isclose(np.array([lon]), np.array([7])).all())

    @pytest.mark.unit
    def test_ease_north_hemi_ease_index2geodetic(self):
        ease = EaseGrid(12000, "NorthHemi")
        x_ind = 500
        y_ind = 400
        lat, lon = ease.ease_index2geodetic(x_ind, y_ind)
        print(lat, lon)
        self.assertTrue(np.isclose(np.array([lat]), np.array([42.419164])).all())
        self.assertTrue(np.isclose(np.array([lon]), np.array([-144.4778])).all())

    @pytest.mark.unit
    def test_ease_north_hemi_geodetic2geodetic(self):
        ease = EaseGrid(12000, "NorthHemi")
        lats = 85
        lons = 7
        (x_ind, y_ind), (x, y) = ease.geodetic2ease(lats, lons)
        lat2, lon2 = ease.ease_coord2geodetic(x, y)
        self.assertTrue(np.isclose(np.array(lats), np.array([lat2])))
        self.assertTrue(np.isclose(np.array(lons), np.array([lon2])))

    @pytest.mark.unit
    def test_ease_global_36km_geodetic2ease(self):
        lons = 17.365144729614258
        lats = 48.57916259765625
        ease_obj = EaseGrid(resolution_m=36000, projection="Global")
        (x_ind, y_ind), (x, y) = ease_obj.geodetic2ease(lat=lats, lon=lons)
        self.assertEqual(528, x_ind)
        self.assertEqual(50, y_ind)

    @pytest.mark.unit
    def test_ease_global_9km_geodetic2ease(self):
        lons = -69.4139
        lats = -22.6355
        ease_obj = EaseGrid(resolution_m=9000, projection="Global")
        (x_ind, y_ind), (x, y) = ease_obj.geodetic2ease(lat=lats, lon=lons)
        self.assertEqual(1184, x_ind)
        self.assertEqual(1124, y_ind)

    @pytest.mark.unit
    def test_ease_north_hemi_36km_geodetic2ease(self):
        lons = 16.9470806121826
        lats = 49.8740196228027
        ease_obj = EaseGrid(resolution_m=36000, projection="NorthHemi")
        (x_ind, y_ind), (x, y) = ease_obj.geodetic2ease(lat=lats, lon=lons)
        self.assertEqual(285, x_ind)
        self.assertEqual(366, y_ind)

    @pytest.mark.unit
    def test_ease_north_hemi_9km_geodetic2ease(self):
        lons = -149.4252
        lats = 69.5271
        ease_obj = EaseGrid(resolution_m=9000, projection="NorthHemi")
        (x_ind, y_ind), (x, y) = ease_obj.geodetic2ease(lat=lats, lon=lons)
        self.assertEqual(871, x_ind)
        self.assertEqual(782, y_ind)

    @pytest.mark.unit
    def test_ease_south_hemi_9km_geodetic2ease(self):
        lons = -69.397
        lats = -22.6699
        ease_obj = EaseGrid(resolution_m=9000, projection="SouthHemi")
        (x_ind, y_ind), (x, y) = ease_obj.geodetic2ease(lat=lats, lon=lons)
        self.assertEqual(264, x_ind)
        self.assertEqual(723, y_ind)

    def tearDown(self):
        pass
