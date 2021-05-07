'''
 ____________________________________________________________________
|TITLE  :  Reading Raster Data with GDAL                             |
|DATE   :  07 May, 2021                                              |
|====================================================================|
|DESCRIPTION:                                                        |
|                                                                    |
|____________________________________________________________________|
'''


#----------------------------#
#== PACKAGES & DIRECTORIES ==#
#----------------------------#
#-- Packages --#
from osgeo import gdal, ogr
import numpy as np
import matplotlib.pyplot as plt
import os

#-- Directories --#
#Set main directory
dirname = os.path.dirname(os.path.realpath("__file__"))

#List paths
paths = dict()
paths["shpfiles_in"] = os.path.join(dirname, "Data", "Inputs", "Shapefiles", "")
paths["rasters_in"] = os.path.join(dirname, "Data", "Inputs",  "Rasters", "")
paths["shpfiles_out"] = os.path.join("Data", "Outputs",  "Shapefiles", "")
paths["rasters_out"] = os.path.join("Data", "Outputs",  "Rasters", "")


os.listdir(os.getcwd())
os.listdir(paths["rasters_in"])

#-----------------------------------#
#== READING & PRE-PROCESSING DATA ==#
#-----------------------------------#
#-- Reading Data --#
ds = gdal.Open(paths["rasters_in"] + "peru_ecoregions.tif")
ds.GetGeoTransform()
