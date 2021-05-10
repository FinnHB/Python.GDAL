'''
 ____________________________________________________________________
|TITLE  :  Reading Raster Data with GDAL                             |
|DATE   :  08 May, 2021                                              |
|====================================================================|
|DESCRIPTION:                                                        |
|                                                                    |
|This script goes through some basic function in GDAL and OGR, using |
|landcover data and a digital elevation model for the United Kingdom.|
|More specifically, the script covers:                               |
|  1. Reading & writing shapefiles/raster files                      |
|  2. Spatial transformations and projections                        |
|  3. Extracting polygons from larger shapefiles                     |
|  4. Clipping raster files based on polygon geometry                |
|  5. Masks raster data by values                                    |
|  6. Converting rasters to polygons                                 |
|  7. Draw new polygons based on coordinates                         |
|____________________________________________________________________|
'''


#----------------------------#
#== PACKAGES & DIRECTORIES ==#
#----------------------------#
#-- Packages --#
import numpy as np
import pandas as pd
from osgeo import gdal, ogr, osr
import matplotlib.pyplot as plt
import os
from itertools import compress



#-- Directories --#
#Set main directory
dirname = os.path.dirname(os.path.abspath("__file__"))
os.listdir(dirname)

#Relative paths
paths = dict()
paths["shp_in"] = os.path.join(dirname, "..", "Data", "Inputs", "Shapefiles", "")
paths["raster_in"] = os.path.join(dirname, "..", "Data", "Inputs",  "Rasters", "")
paths["shp_out"] = os.path.join(dirname, "..", "Data", "Outputs",  "Shapefiles", "")
paths["raster_out"] = os.path.join(dirname, "..", "Data", "Outputs",  "Rasters", "")



#----------------#
#== PARAMETERS ==#
#----------------#
#-- Switches --#
overwrite = False                                                     #If False, will not overwrite data if it already exists


#-----------------------------------#
#== READING & PRE-PROCESSING DATA ==#
#-----------------------------------#
#-- Reading Data --#
#Reading corine land cover
corine = gdal.Open(paths["raster_in"] + "CLC_2018.tif")

#UK border & counties
uk = ogr.Open(paths["shp_in"] + "uk_poly.shp", 1)

#Combine all the DEMs for North Yorkshire into a single DEM
ny_dem = gdal.Open(paths["raster_in"] + "nyorkshire_dem.tif")


#-- Reproject to WGS84 --#
#Corine land cover
fname = "corine_wgs84.tif"
if fname not in os.listdir(paths["raster_in"]) or overwrite == True:
    file = paths["raster_out"] + fname
    corine_wgs84 = gdal.Warp(file, corine, dstSRS = "EPSG:4326")

#North Yorkshire DEM
fname = "nyorkshire_dem_wgs84.tif"
if fname not in os.listdir(paths["raster_in"]) or overwrite == True:
    file = paths["raster_out"] + fname
    dem_wgs84 = gdal.Warp(file, ny_dem, dstSRS = "EPSG:4326")


#-- Reduce resolution --#
#Corine
file = paths["raster_out"] + "corine_lowres.tif"
infile = paths["raster_out"] + "corine_wgs84.tif"
corine_lr = gdal.Warp(file, infile, xRes = 0.01, yRes = 0.01,
                      resampleAlg = "min")
#Clip landcover to UK
fname = "corine_uk.tif"
if fname not in os.listdir(paths["raster_out"]) or overwrite == True:
    #Reads in Corine data if not already loaded
    file = paths["raster_out"] + fname
    crop_poly = paths["shp_in"] + "uk_poly.shp"
    corine_uk = gdal.Warp(file, corine_lr, cutlineDSName = crop_poly,
    cropToCutline = True, dstNodata = np.nan)

#DEM
file = paths["raster_out"] + "dem_lowres.tif"
infile = paths["raster_out"] + "nyorkshire_dem_wgs84.tif"
dem_lr = gdal.Warp(file, infile, xRes = 0.01, yRes = 0.01,
                      resampleAlg = "min")
