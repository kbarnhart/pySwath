# -*- coding: utf-8 -*-
"""
Created on Fri Feb 07 16:53:28 2014

@author: Lab User
"""

# Name: FocalStatistics_Ex_02.py
# Description: Calculates a statistic on a raster over a specified
#    neighborhood.
# Requirements: Spatial Analyst Extension
# Author: ESRI

# Import system modules
import arcpy
from arcpy import env
from arcpy.sa import *

# Set environment settings
env.workspace = r"C:\Users\Lab User\Billy\swath_profiler"

# Set local variables
inRaster = "julJul2013_corr3_dailyDisp_forArc.tif"
neighborhood = NbrRectangle(4, 4, "CELL")

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

# Execute FocalStatistics
outFocalStatistics = FocalStatistics(inRaster, neighborhood, "MEAN","DATA")

# Save the output 
outFocalStatistics.save("C:\Users\Lab User\Billy\swath_profiler\name")