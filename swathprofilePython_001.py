### Swath profiler code ###
# Written by William Armstrong and Katy Barnhart, with code modified from KRB's Arc-based profiler
# Begun 15 July 2015
# Finished xxxx

print "STARTING SWATH PROFILER"
print "Importing modules"    

# ultimately we want to make this a function, with the following input information:

##### USER DEFINED VARIABLES ######

folderPath = u'/Users/wiar9509/git/pySwath'
inPoly='/test/testMultiLine.shp' # input shapefile filename along which raster will be sampled
inRast='/test/vv.tif' # input raster filename to sample
width=500.0 # across-flow box dimension [this is in the units of the projection]
height=150.0 # along-flow box dimension

outTextFilename='test.csv' # filename for textfile output
fOutFigure='fOut.pdf' # filename for figure output

###### LOADING MODULES ######

import os
os.chdir(folderPath)
import numpy as np
from krb_vecTools import *
#import matplotlib.pyplot as plt
import scipy as sci 
import grass.script as grass

from osgeo import ogr # modify code to just use ogr, not arcGIS
from osgeo import osr


###### PROCESSING ######

tempPoly='temp.shp'
tempPoly2='temp2.shp'
tempPoly3='temp3.shp'

# create the spatial reference
srs = osr.SpatialReference()
srs.ImportFromEPSG(32607) # This hard-coded to UTM 7N w/ WGS84 datum. Would be nice to have this defined based of input shapefile

print "Opening shapefile with OGR"
shapefile=folderPath+inPoly
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(shapefile, 0)
layer = dataSource.GetLayer()
crs = layer.GetSpatialRef() # Coordinate reference system

numLines=layer.GetFeatureCount() # Number of lines in shapefile

# Iterating over number of lines in input shapefile
for lineNo in range(numLines):
	
	# Initializing
	lineEast=[]
	lineNorth=[]
	dE=[]
 	dN=[]
 	sampleEast=[]
 	sampleNorth=[]
 	ptArray=[]
 	
 	feat=layer.GetFeature(lineNo) # Highlights current line
	geom=feat.geometry()
	name=feat.GetField(1) # Gets current line name
	
	print "PROCESSING TRANSECT "+name.upper()
	print "TRANSECT "+str(lineNo+1)+" OF "+str(numLines)
	print "Getting coordinates along line"
	
	numPoints=geom.GetPointCount() # Number of vertices in polyline
	
	for nP in range(numPoints): # Iterates over points and appends x,y (easting, northing) coords from current line
		lineEast.append(geom.GetPoint(nP)[0]) # change in easting from last coordinate [m]
		lineNorth.append(geom.GetPoint(nP)[1]) # change in northing from last coordinate [m]
		if nP > 0: # Calculates difference in easting/northing between points for later line length calculation
			dE.append(lineEast[nP]-lineEast[nP-1])
			dN.append(lineNorth[nP]-lineNorth[nP-1])
	
	# Converting to arrays to do math
	dE=np.array(dE)
	dN=np.array(dN)
	
	# Finding centroid coordinates for sampling
	lineSegs=len(dE) # number of line segments
	segmentLen=np.sqrt(dE**2 + dN**2) # units same as projection
	pointsInSeg=np.ceil(segmentLen/height)
	dEseg=dE/pointsInSeg
	dNseg=dN/pointsInSeg
	
	# Iterate over line segments to generate coordinates
	for segNo in range(lineSegs):
		index=np.linspace(1,pointsInSeg[segNo]+1,pointsInSeg[segNo]+1) # list from 1 to number of sampling points in line segment
		sampleEast.append(lineEast[segNo]+dEseg[segNo]*index)
		sampleNorth.append(lineNorth[segNo]+dNseg[segNo]*index)
		#plt.plot(sampleEast[segNo],sampleNorth[segNo],'.')
		
	# Was in a weird format before because of segNo loop. This makes it a better array
	sampleEast=np.concatenate(sampleEast,axis=0)
	sampleNorth=np.concatenate(sampleNorth,axis=0)
	
	# Converting from array to list
	sampleEast=sampleEast.tolist()
	sampleNorth=sampleNorth.tolist()
	
	# Putting easting and northing coordinates together into a list of lists
	for coord in range(len(sampleEast)):
		ptArray.append((sampleEast[coord],sampleNorth[coord]))
	
	# Plotting to make sure everything working	
	#plt.plot(lineEast,lineNorth,'o')
	#plt.show()
	#plt.axis('equal')
	
	# Vectools to generate polygons (rectangles) for zonal stats 
	print "Calculating the slope and perpendicular slope along the line"
	hwin=5
	slope, perpSlope=calcSlopeAndPerpSlope(ptArray, hwin)
	#slope, perpSlope=calcSlopeAndPerpSlope(ptArray,int(height))

	print "Storing the coordinates of polygons for sampling"
	polygons=[]
	for i in range(len(ptArray)):
		polygons.append(makePoly(ptArray[i], perpSlope[i], width, height))
	# Format of polygon coordinates is lower right (LR), LL, UL, UR

	#calculate distance along the line
 	lineDist=distanceAlongLine(sampleEast,sampleNorth)
	
	### CREATING POLYGONS ###
	
	# Sourced ideas for the below lines from http://www.gis.usu.edu/~chrisg/python/2008/os2_slides.pdf as well as https://pcjericks.github.io/py-gdalogr-cookbook/

	# Initializing
	numBoxes=len(polygons)

	print "Generating and populating polygons"
	# I am not sure what a lot of this means, but it is required to build a shapefile
	if os.path.exists(name+"polygons.shp"):
		driver.DeleteDataSource(name+"polygons.shp") # error if data source already exists
	newDataSource=driver.CreateDataSource(name+"polygons.shp")
	newLayer=newDataSource.CreateLayer(name+"polygons",srs,geom_type=ogr.wkbPolygon)
	fieldDefn=ogr.FieldDefn('id',ogr.OFTInteger)
	#newLayer.CreateField(ogr.FieldDefn("name",ogr.OFTString))
	newLayer.CreateField(fieldDefn)
	newLayer.CreateField(ogr.FieldDefn("east",ogr.OFTReal))
	newLayer.CreateField(ogr.FieldDefn("north",ogr.OFTReal))
	newLayer.CreateField(ogr.FieldDefn("dist",ogr.OFTReal))
	newLayer.CreateField(ogr.FieldDefn("min",ogr.OFTReal))
	newLayer.CreateField(ogr.FieldDefn("mean",ogr.OFTReal))
	newLayer.CreateField(ogr.FieldDefn("max",ogr.OFTReal))
	newLayer.CreateField(ogr.FieldDefn("range",ogr.OFTReal))
	newLayer.CreateField(ogr.FieldDefn("stddev",ogr.OFTReal))
	newLayer.CreateField(ogr.FieldDefn("sum",ogr.OFTReal))

	

	# Iterating to create polygons and set fields
	for poly in range(numBoxes):
		# create the feature
		feature = ogr.Feature(newLayer.GetLayerDefn())
		
		# Set Geometry
		ring=ogr.Geometry(ogr.wkbLinearRing)
		ring.AddPoint(polygons[poly][0][0],polygons[poly][0][1]) # adding easting/northing for each vertex
		ring.AddPoint(polygons[poly][1][0],polygons[poly][1][1])
		ring.AddPoint(polygons[poly][2][0],polygons[poly][2][1])
		ring.AddPoint(polygons[poly][3][0],polygons[poly][3][1])
		ring.CloseRings()

		polygon = ogr.Geometry(ogr.wkbPolygon)
		polygon.AddGeometry(ring)
		feature.SetGeometryDirectly(polygon)
		
		# Setting fields
		feature.SetField('id',poly) # id number
		featPoly=polygon.GetGeometryRef(0) # pointer
		featEast=featPoly.GetX()
		featNorth=featPoly.GetY()
		feature.SetField('east',featEast) # easting (centroid I think)
		feature.SetField('north',featNorth) # northing
		featDist=lineDist[poly]
		feature.SetField('dist',featDist) # distance along line
		
		# Create the feature in the layer (shapefile)
		newLayer.CreateFeature(feature)
		
		# Destroy the feature to free resources
		feature.Destroy()

	# Destroy the data source to free resources
	newDataSource.Destroy()

	#### GRASS PORTION ####
	print "Creating features in GRASS and performing zonal statistics"

	# Read in raster to sample
	grass.run_command("r.in.gdal",flags='e',overwrite=True,input=folderPath+inRast,output='sampleRast')
	# Read in polygons to use as sampling bins
	grass.run_command("v.in.ogr",overwrite=True,input=folderPath+'/'+name+"polygons.shp",output='samplePolys')
	# Make sure computational region contains polygons
	grass.run_command("g.region",vec='samplePolys',res='10')
	# Convert polygons vector into 'zone' raster for sampling
	grass.run_command("v.to.rast",overwrite=True,input='samplePolys',output='zoneRast',use='attr',attr='id')
	# Run univariate statistics on sample raster, using zone raster to bin and populate new rows
	grass.run_command("r.univar",overwrite=True,flags='t',map='sampleRast', zones='zoneRast', out=outTextFilename, sep=',',)
	# You now have a file named $outTextFilename in the folder in which this swath profiler code lives
	
	### END GRASS PORTION ###
	
	min=[]
	max=[]
	rangeVals=[]
	mean=[]
	stdDev=[]
	sum=[]
	numPixels=[]
	
	# Read in stats file
	print 'Reading in statistics file and updating shapefile'
	statsFile=np.genfromtxt(outTextFilename,delimiter=',')
	statsLen=len(statsFile)
	
	# Opening layer to populate stats values in shapefile created above
	updateShapefile=ogr.Open(name+"polygons.shp")		
	updateLayer=updateShapefile.GetLayer(0)
	
	for line in range(statsLen): # First value of these is nan
		min.append(statsFile[line][4]) # minimum value of raster within polygon
		max.append(statsFile[line][5]) # maximum value
		rangeVals.append(statsFile[line][6]) # range of values
		mean.append(statsFile[line][7]) # mean value
		stdDev.append(statsFile[line][9]) # standard deviation of values
		sum.append(statsFile[line][12])	# sum of values

	
	# 	print 'Save Values'
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
