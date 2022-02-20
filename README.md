# easepy

![test-main](https://github.com/karl-nordstrom/easepy/actions/workflows/python-test-main.yml/badge.svg)
![coverage-main](https://img.shields.io/codecov/c/github/karl-nordstrom/easepy)
![license](https://img.shields.io/github/license/karl-nordstrom/easepy)

A python package for working with EASE grids in geodetic coordinates.
The documentation is available at https://easepy.readthedocs.io/en/latest/.

Installation
------------

    pip install easepy

Example usage
-------------

    import easepy
    ease = easepy.EaseGrid(resolution_m=25000, projection="Global")
    # Fetch grid cell centroids
    grid_lats, grid_lons = ease.geodetic_grid
    # Find corresponding cell indices for particular location(s)
    ease_indices, _ = ease.geodetic2ease(lat=46.7, lon=132.1)

Authors:

- Karl Nordstrom (<karl.am.nordstrom@gmail.com>)
- Giorgio Savastano (<giorgiosavastano@gmail.com>)

Please use github issues to make bug reports and request new functionality. Contributions are always welcome.
