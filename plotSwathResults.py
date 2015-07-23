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
	
#	print 'Save Values'
# 	saveOut=[np.asarray(xOut), np.asarray(yOut), np.asarray(distOut),np.asarray(meanOut), np.asarray(minOut), np.asarray(maxOut)] 
# 	np.savetxt(fOutText, np.asarray(saveOut).T, delimiter=',')
# 
# 	print "Make Plot"
# 	fig1=plt.figure(num=None, figsize=(3.35, 5), dpi=300, facecolor='w', edgecolor='w')
# 	#fig1.patch.set_alpha(0.0)
# 
# 	ax = fig1.add_subplot(1,1,1)
# 
# 	p1,=plt.plot(distOut, meanOut, lw=2, color='k', ls='-')
# 	p2,=plt.plot(distOut, minOut, lw=1, color='grey', ls='-')
# 	p3,=plt.plot(distOut, maxOut, lw=1, color='grey', ls='--')
# 
# 	lg1=plt.legend([p1, p2, p3], ['Mean', 'Minimum', 'Maximum'], loc='upper left', fancybox=True,bbox_to_anchor=(0, 1))
# 	lg1.get_frame().set_edgecolor('grey')
# 	lg1.get_frame().set_linewidth(0.5)
# 
# 	plt.tick_params(axis='both', which='major', labelsize=6)
# 	plt.tick_params(axis='both', which='minor', labelsize=6)
# 
# 	plt.xlabel('Distance along Profile',fontsize=8)    
# 	plt.ylabel('Value',fontsize=8) 
# 	plt.title('Swath' ,fontsize=10)
# 
# 	plt.savefig(fOutFigure, format='pdf')
# 	plt.show()
