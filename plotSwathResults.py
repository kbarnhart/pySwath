### SCRIPT TO READ IN AND PLOT SWATH PROFILER RESULTS
# Written by William Armstrong
# 23 July 2015

### USER DEFINED VARIABLES ###
folderPath = u'/Users/wiar9509/git/pySwath'
inPolySource='/test/testMultiLine.shp' # initial shapefile filename used for swath profiler
inRast='test.tif' # input raster filename from swath profiler
#inPoly= 0 # to read a specific shapefile
###### LOADING MODULES ######

import os
os.chdir(folderPath)
import numpy as np
import matplotlib.pyplot as plt

from osgeo import ogr # modify code to just use ogr, not arcGIS
from osgeo import osr
	
###### PROCESSING ######

# Function to plot results given input transect name
def plotSwath(polyName):
	print 'test'
	

# Read in shapefiles
print "Opening shapefile with OGR"
shapefile=folderPath+inPolySource
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(shapefile, 0)
layer = dataSource.GetLayer()
crs = layer.GetSpatialRef() # Coordinate reference system

numLines=layer.GetFeatureCount() # Number of lines in shapefile

for line in range(numLines):

	# Initializing
	east=[]
	north=[]
	dist=[]
	min=[]
	mean=[]
	max=[]
	stdDev=[]
	
	trans=layer.GetFeature(line)
	transectName=trans.GetField(1)
	shapeName=transectName+'polygons.shp'
	
	# Opening layer to plot
	shapeDatasource = driver.Open(shapeName, update=0)
	shapeLayer = shapeDatasource.GetLayer()
	
	numBoxes = shapeLayer.GetFeatureCount()
	
	for box in range(numBoxes):
		boxNow=shapeLayer.GetFeature(box)
		if boxNow.GetField(5) is not None: # some are 'None'
			east.append(boxNow.GetField(1))
			north.append(boxNow.GetField(2))
			dist.append(boxNow.GetField(3))
			min.append(boxNow.GetField(4))
			mean.append(boxNow.GetField(5))
			max.append(boxNow.GetField(6))
			stdDev.append(boxNow.GetField(8))
			
	dist=np.array(dist)
	
	titleText=transectName+' profile from raster: '+inRast
	ax=plt.subplot(111)
	ax.plot(dist/1000,mean,color='k',linewidth=2,label='Mean')
	ax.fill_between(dist/1000.0,min,max,color='gray',alpha=0.5,label='Range')
	ax.set_title(titleText,fontsize=18)
	ax.set_xlabel('Distance [km]',fontsize=14)
	ax.set_ylabel('Velocity [m d$^{-1}$]',fontsize=14)
	ax.legend(loc='best')
	ax.grid()
	plt.savefig(inRast+'_profile.pdf',orientation='landscape',format='pdf')
	plt.show()
	plt.draw()