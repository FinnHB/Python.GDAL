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
import subprocess



#-- Directories --#
#Set main directory
dirname = os.path.dirname(os.path.abspath("__file__"))
os.listdir(dirname)

#Relative paths
paths = dict()
paths["shp_in"] = os.path.join(dirname, "..", "Data", "Inputs", "Shapefiles", "")
paths["raster_in"] = os.path.join(dirname, "..", "Data", "Inputs",  "Rasters", "")
paths["data_in"] = os.path.join(dirname, "..", "Data", "Inputs",  "Datasets", "")
paths["shp_out"] = os.path.join(dirname, "..", "Data", "Outputs",  "Shapefiles", "")
paths["raster_out"] = os.path.join(dirname, "..", "Data", "Outputs",  "Rasters", "")
paths["images"] = os.path.join(dirname, "..", "Images", "")

#-- Notes --#
# Corine land-cover classifications can be found here:
# https://land.copernicus.eu/user-corner/technical-library/corine-land-cover-nomenclature-guidelines/html



#----------------#
#== PARAMETERS ==#
#----------------#
#-- Switches --#
overwrite = False                                                     #If False, will not overwrite data if it already exists
show_plots = True                                                     #If True, shows plots
save_plots = True                                                     #If True, saves DEM plot of North Yorkshire (cannot save without showing plots)

#-- Analysis Parameters --#
#DEM Data Settings
elevation_region = "mid"                                              #Specify "low", "mid" or "high" elevation areas (will extracting raster values)

#Corine Data Settings
corine_regions = (211,245)                                            #Lower and uppbound land cover classification (211-244 = Agriculture)
corine_nomclature = "agriculture"                                     #Name to be used to identify tif file subset by corine_regions (not including file suffix)



#-----------------------------------#
#== READING & PRE-PROCESSING DATA ==#
#-----------------------------------#
#-- Reading Data --#
#Reading corine land cover
corine = gdal.Open(paths["raster_in"] + "corine_lc_2018.tif")

#UK border & counties
uk = ogr.Open(paths["shp_in"] + "uk_poly.shp", 1)                     #0 = immutable, 1 = mutable
counties = ogr.Open(paths["shp_in"] + "uk_counties.shp", 1)

#Combine all the DEMs for North Yorkshire into a single DEM
ny_dem = gdal.Open(paths["raster_in"] + "nyorkshire_dem.tif")


#-- Reproject to WGS84 --#
#Corine land cover
fname = "corine_wgs84.tif"
if fname not in os.listdir(paths["raster_out"]) or overwrite == True:
    file = paths["raster_out"] + fname
    corine_wgs84 = gdal.Warp(file, corine, dstSRS = "EPSG:4326")

#North Yorkshire DEM
fname = "nyorkshire_dem_wgs84.tif"
if fname not in os.listdir(paths["raster_out"]) or overwrite == True:
    file = paths["raster_out"] + fname
    dem_wgs84 = gdal.Warp(file, ny_dem, dstSRS = "EPSG:4326")



#-- Extract County of North Yorkshire --#
#Select polygon
counties_layer = counties.GetLayer()
ny_oid = np.where([f.GetField("ctyua16nm") == "North Yorkshire" for f in counties_layer])
ny_oid = int(ny_oid[0][0])
nyfeature = counties_layer.GetFeature(ny_oid)
nyname = nyfeature.GetField("ctyua16nm")
nygeom = nyfeature.GetGeometryRef()

#Save shapefile
outfile = paths["shp_out"] + "nyorkshire.shp"
driver = ogr.GetDriverByName("ESRI Shapefile")
outds = driver.CreateDataSource(outfile)
layer = outds.CreateLayer("north_yorkshire", counties_layer.GetSpatialRef())
feature = ogr.Feature(layer.GetLayerDefn())
feature.SetGeometry(nygeom)
layer.CreateFeature(feature)

#Destroy/close shapefile, this prevents the file from getting corrupted.
nygeom.Destroy()
outds.Destroy()



#-- Clipping Data to North Yorkshire --#
#Clipping corine Corine Landcover
fname = "corine_ny.tif"
if fname not in os.listdir(paths["raster_out"]) or overwrite == True:
    #Reads in Corine data if not already loaded
    if "corine_wgs84" not in locals():
        corine_wgs84 = gdal.Open(paths["raster_out"] + "corine_wgs84.tif")
    file = paths["raster_out"] + fname
    crop_poly = paths["shp_out"] + "nyorkshire.shp"
    corine_ny = gdal.Warp(file, corine_wgs84, cutlineDSName = crop_poly,
                          cropToCutline = True, dstNodata = np.nan)

#Clipping DEM
fname = "dem_ny.tif"
if fname not in os.listdir(paths["raster_out"]) or overwrite == True:
    #Reads in DEM if not already loaded
    if "corine_wgs84" not in locals():
        dem_wgs84 = gdal.Open(paths["raster_out"] + "nyorkshire_dem_wgs84.tif")
    file = paths["raster_out"] + fname
    crop_poly = paths["shp_out"] + "nyorkshire.shp"
    dem_ny = gdal.Warp(file, dem_wgs84, cutlineDSName = crop_poly,
                       cropToCutline = True, dstNodata = np.nan)



#-- Visualising Results --#
#Loading data if not loaded
if "corine_ny" not in locals():
    corine_ny = gdal.Open(paths["raster_out"] + "corine_ny.tif")
if "dem_ny" not in locals():
    dem_ny = gdal.Open(paths["raster_out"] + "dem_ny.tif")

#Getting values
corine_vals = corine_ny.GetRasterBand(1).ReadAsArray()
dem_vals = dem_ny.GetRasterBand(1).ReadAsArray()

if show_plots:
    #Image of both maps
    fig, ax = plt.subplots(1, 2, figsize=(24, 12))
    ax[0].title.set_text("CORINE Land Cover")
    ax[0].imshow(corine_vals)
    ax[1].title.set_text("Digital Elevation Model")
    ax[1].imshow(dem_vals)

    #Digital Elevation Map of North Yorkshire
    if save_plots:
        plt.figure(figsize = [18,18])
        plt.imshow(dem_vals, cmap = "RdBu_r")
        plt.title("Digital Elevation Map of North Yorkshire", fontsize = 25)
        cbar = plt.colorbar(shrink = 0.3)
        cbar.set_label("Elevation (m)")
        plt.savefig(paths["images"] + "north_yorkshire_dem.png")





#---------------------------#
#== MASK BY RASTER VALUES ==#
#---------------------------#
#-- Elevation --#
#Identify elevation level based on parameters set at start of script
dem_mean = np.nanmean(dem_vals)
dem_sd = np.nanstd(dem_vals)
if elevation_region == "low":
    mask = np.where((dem_vals < (dem_mean-dem_sd)),1,0)
elif elevation_region == "high":
    mask = np.where((dem_vals > (dem_mean+dem_sd)),1,0)
else:
    upper_bound = dem_vals < (dem_mean+dem_sd)
    lower_bound = dem_vals > (dem_mean-dem_sd)
    mask = np.where((lower_bound & upper_bound),1,0)


#Saving the mask as a raster
fname = "ny_{}lands.tif".format(elevation_region)
if fname not in os.listdir(paths["raster_out"]) or overwrite == True:
    file = paths["raster_out"] + fname
    driver = gdal.GetDriverByName("GTiff")
    driver.Register()
    outds = driver.Create(file, xsize = mask.shape[1],
                          ysize = mask.shape[0], bands = 1,
                          eType = gdal.GDT_Int16)
    outds.SetGeoTransform(dem_ny.GetGeoTransform())
    outds.SetProjection(dem_ny.GetProjection())
    outband = outds.GetRasterBand(1)
    outband.WriteArray(mask)
    outband.SetNoDataValue(np.nan)
    outband.FlushCache()
    outband = outds = None



#-- Agricultural areas --#
#Identify mask based on corine regions
upper_bound = corine_vals < max(corine_regions)
lower_bound = corine_vals > min(corine_regions)
masked_area = upper_bound & lower_bound
mask = np.where(masked_area,1,0)


#Saving the mask as a raster
fname = corine_nomclature + ".tiff"
if fname not in os.listdir(paths["raster_out"]) or overwrite == True:
    file = paths["raster_out"] + fname
    driver = gdal.GetDriverByName("GTiff")
    driver.Register()
    outds = driver.Create(file, xsize = mask.shape[1],
                          ysize = mask.shape[0], bands = 1,
                          eType = gdal.GDT_Int16)
    outds.SetGeoTransform(corine_ny.GetGeoTransform())
    outds.SetProjection(corine_ny.GetProjection())
    outband = outds.GetRasterBand(1)
    outband.WriteArray(mask)
    outband.SetNoDataValue(np.nan)
    outband.FlushCache()
    outband = outds = None



#-- Reading in the files just made --#
#Corine land use
file = paths["raster_out"] + "ny_{}lands.tif".format(elevation_region)
dem = gdal.Open(file)

#DEM
file = paths["raster_out"] + corine_nomclature + ".tiff"
corine_lu = gdal.Open(file)


#Create version with the same resolution as the corine data
gt = corine_lu.GetGeoTransform()
gt = dem.GetGeoTransform()



#-- Visualising Results --#
#Getting values
corine_vals = corine_lu.GetRasterBand(1).ReadAsArray()
dem_vals = dem.GetRasterBand(1).ReadAsArray()

#Image of the maps
if show_plots:
    fig, ax = plt.subplots(1, 2, figsize=(24, 14))
    ax[0].title.set_text("CORINE Land Cover: {}".format(corine_nomclature))
    ax[0].imshow(corine_vals)
    ax[1].title.set_text("Digital Elevation Model: binary {}elevation".format(elevation_region))
    ax[1].imshow(dem_vals)
    if save_plots:
        plt.savefig(paths["images"] + "binary_mask.png")







#-------------------#
#== DRAW POLYGONS ==#
#-------------------#
#-- Polygonises Rasters --#
fname = "dem_polygon.shp"
if fname not in os.listdir(paths["raster_out"]) or overwrite == True:
    file = paths["shp_out"] + fname
    dem_band = dem.GetRasterBand(1)
    driver = ogr.GetDriverByName("ESRI Shapefile")
    outds = driver.CreateDataSource(file)
    layer = outds.CreateLayer("polygonized", srs=None)
    gdal.Polygonize(dem_band, None, layer, -1, [], callback=None)
    outds.Destroy()
    sourceRaster = None


#-- Draw polygon from coordinates --#
#Get extent of corine raster
gt = dem.GetGeoTransform()
x_min = gt[0]
y_max = gt[3]
x_max = x_min + gt[1] * dem.RasterXSize
y_min = y_max + gt[5] * dem.RasterYSize

#Create a boundary box
square = ogr.Geometry(ogr.wkbLinearRing)
square.AddPoint(x_min, y_max)
square.AddPoint(x_max, y_max)
square.AddPoint(x_max, y_min)
square.AddPoint(x_min, y_min)
square.AddPoint(x_min, y_max)

#Write boundary box
filename = paths["shp_out"] + "simple_box.shp"
poly = ogr.Geometry(ogr.wkbPolygon)
poly.AddGeometry(square)
driver = ogr.GetDriverByName("ESRI Shapefile")
outds = driver.CreateDataSource(filename)
outlayer = outds.CreateLayer("box", dem.GetSpatialRef())
feature = ogr.Feature(outlayer.GetLayerDefn())
feature.SetGeometry(poly)
outlayer.CreateFeature(feature)
outds = outlayer = feature = None


#-- Draw polygons from csv --#
#Read csv
file = paths["data_in"] + "poly_coords.csv"
coords = pd.read_csv(file)
poly_ids = np.unique(coords.pid)

#Initialise polygon
filename = paths["shp_out"] + "unsquare_box.shp"
poly = ogr.Geometry(ogr.wkbPolygon)

#Looping over polygon id's
for pid in poly_ids:
    #Subset coords
    p_coords = coords.loc[coords.pid==pid,:]

    #Creating polygon
    square = ogr.Geometry(ogr.wkbLinearRing)
    [square.AddPoint(x,y) for x,y in p_coords[["x","y"]].values]
    square.CloseRings()                         #Ensuring that the polygon ends where it started

    #Assign to polygon
    poly.AddGeometry(square)

#Write polygon
driver = ogr.GetDriverByName("ESRI Shapefile")
outds = driver.CreateDataSource(filename)
outlayer = outds.CreateLayer("box", dem.GetSpatialRef())
feature = ogr.Feature(outlayer.GetLayerDefn())
feature.SetGeometry(poly)
outlayer.CreateFeature(feature)
outds = outlayer = feature = None
