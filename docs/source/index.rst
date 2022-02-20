.. easepy documentation master file, created by
   sphinx-quickstart on Sat Feb 19 15:44:06 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to easepy's documentation!
==================================

easepy allows for gridding into the Equal Area Scaleable Earth (EASE) grids.
These grids are popular for global-scale earth observation data. More details
on the EASE grids can be found at https://nsidc.org/data/ease.

To install:

.. code-block:: console

   pip install easepy

Example usage
-------------

.. code-block:: python

   import easepy
   ease = easepy.EaseGrid(resolution_m=25000, projection="Global")
   # Fetch grid cell centroids
   grid_lats, grid_lons = ease.geodetic_grid
   # Find corresponding cell indices for particular location(s)
   ease_indices, _ = ease.geodetic2ease(lat=46.7, lon=132.1)

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. autoclass:: easepy.EaseGrid
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
