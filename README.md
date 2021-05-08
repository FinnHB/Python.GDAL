Python.GDAL
## GDAL with the Python API ##
#### Core functions to get starter working with spatial data ####

This repository was set up in particular to explore the functionalities of GDAL and OGR. The directory contains a script which looks at the reads in several raster and polygon layers, cropping them to the extent of North Yorkshire. As of yet, the repository contains a single script which covers some of the basic functions of GDAL, namely:

1. Reading and writing shapefiles and rasters
2. Spatial transformations and projections
3. Extracting features from larger polygons
4. Clipping raster files by polygon geometry
5. Masking raster data by values
6. Converting raster data to polygons
7. Creating new polygons based on lists of coordinates.

The resulting map of the cropped Digital Elevation Model can be seen below.

![plot](https://github.com/FinnHB/Python.GDAL/blob/main/Images/north_yorkshire_dem.png)

A binary (0,1) mask was also constructed based on the a given elevation range and land-use type.

![plot](https://github.com/FinnHB/Python.GDAL/blob/main/Images/binary_mask.png)
