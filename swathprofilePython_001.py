print "importing modules"    


# ultimately we want to make this a function, with the following input information:

### USER DEFINED VARIABLES ####

folderPath = u'/Users/wiar9509/git/pySwath'
inPoly='/test/testMultiLine.shp' # input shapefile filename along which raster will be sampled
inRast='vv.tif' # input raster filename to sample
width=100.0 # across-flow box dimension [this is in the units of the projection]
height=50.0 # along-flow box dimension

fOutText='fOutCSV.csv' # filename for textfile output
fOutFigure='fOut.pdf' # filename for figure output

### LOADING MODULES ###

import os
os.chdir(folderPath)
import numpy as np
from krb_vecTools import *
import matplotlib.pyplot as plt
import scipy as sci 

from osgeo import ogr # modify code to just use ogr, not arcGIS
from osgeo import osr


### PROCESSING ###

tempPoly='temp.shp'
tempPoly2='temp2.shp'
tempPoly3='temp3.shp'

print "Opening shapefile with OGR"
shapefile=folderPath+inPoly
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(shapefile, 0)
layer = dataSource.GetLayer()
crs = layer.GetSpatialRef() # Coordinate reference system

numLines=layer.GetFeatureCount() # Number of lines in shapefile

for lineNo in range(numLines): # need to add additional info for output files if there are more than one line. 20jul15 WHA - did this in 'name' variable below
	
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
	
	print "Getting coordinates along line: "+name
	
	numPoints=geom.GetPointCount()
	
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
	
	for segNo in range(lineSegs): # iterate over line segments to generate coordinates
		index=np.linspace(1,pointsInSeg[segNo]+1,pointsInSeg[segNo]+1) # list from 1 to number of sampling points in line segment
		sampleEast.append(lineEast[segNo]+dEseg[segNo]*index)
		sampleNorth.append(lineNorth[segNo]+dNseg[segNo]*index)
		plt.plot(sampleEast[segNo],sampleNorth[segNo],'.')
		
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
	plt.plot(lineEast,lineNorth,'o')
	plt.show()
	plt.axis('equal')
	
	# Vectools to generate polygons (rectangles) for zonal stats 
	print "Calculating the slope and perpendicular slope along the line"
	hwin=5
	slope, perpSlope=calcSlopeAndPerpSlope(ptArray, hwin)

	print "Storing the coordiates of the swath cross section lines"
	polygons=[]
	for i in range(len(ptArray)):
		polygons.append(makePoly(ptArray[i], perpSlope[i], width, height))
	# Format of polygon coordinates is lower right (LR), LL, UL, UR
	
	### CREATING POLYGONS ###
	
	# Sourced ideas for the below lines from http://www.gis.usu.edu/~chrisg/python/2008/os2_slides.pdf
	
	# Initializing
	multiBox= ogr.Geometry(ogr.wkbMultiPolygon)

	# Iterating over coordinates to create polygons
	for poly in range(len(polygons)):
		ring=ogr.Geometry(ogr.wkbLinearRing)
		ring.AddPoint(polygons[poly][0][0],polygons[poly][0][1]) # adding easting/northing for each vertex
		ring.AddPoint(polygons[poly][1][0],polygons[poly][1][1])
		ring.AddPoint(polygons[poly][2][0],polygons[poly][2][1])
		ring.AddPoint(polygons[poly][3][0],polygons[poly][3][1])
		ring.CloseRings()
		
		box=ogr.Geometry(ogr.wkbPolygon)
		box.AddGeometry(ring)
		multiBox.AddGeometry(box)
	
	# I am not sure what a lot of this means, but it is required to build a shapefile
	fieldDefn=ogr.FieldDefn('id',ogr.OFTInteger)
	if os.path.exists(name+"multiBox.shp"):
		driver.DeleteDataSource(name+"multiBox.shp") # error if data source already exists
	newDataSource=driver.CreateDataSource(name+"multiBox.shp")
	newLayer=newDataSource.CreateLayer('test',geom_type=ogr.wkbMultiPolygon)
	newLayer.CreateField(fieldDefn)
	featureDefn=newLayer.GetLayerDefn()
	feature=ogr.Feature(featureDefn)
	feature.SetGeometry(multiBox) # not sure if this right; outputs '0'
	feature.SetField('id',1)
	newLayer.CreateFeature(feature)
	
	newDataSource.Destroy() # closes data source and writes shapefile
	
	
	#fieldDefn.SetWidth(5)

# 		polyNow = 'box'+str(poly)
# 		eval(polyNow+" = ogr.Geometry(ogr.wkbLinearRing)")
# 		eval(polyNow+".addPoint(polygons[0])")
# 		polynow.addPoint(polygons[1])
# 		polynow.addPoint(polygons[2])
# 		polynow.addPoint(polygons[3])

		#box+poly = (1,1)
# 	#calculate distance along the line
# 	distOut=distanceAlongLine(sampleEast,sampleNorth)
# 
# 	# initialize spatial containers (Need to revise with OGR output)
# 	# Below lines copied from "Create a new shapefile and add data" on https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#create-a-new-shapefile-and-add-data. I am pretty confused about what each step was actually doing.
# 	
# 	# Create data source
# 	newShapefile=driver.CreateDataSource(name+"_polygons.shp")
# 	
# 	# Create spatial reference
# 	srs=osr.SpatialReference()
# 	srs.ImportFromEPSG(32607) # This hard-coded to UTM 7N w/ WGS84 datum. Would be nice to have this defined based of input shapefile
# 	
# 	# Create layer
# 	newLayer=newShapefile.CreateLayer(name,srs,ogr.wkbPoint)
# 	
# 	# Adding fields to shapefile
# 	boxNum=ogr.FieldDefn("boxNum",ogr.OFTInteger)
# 	boxNum
# 	
	
	# point=arcpy.Point()
# 	array=arcpy.Array()
# 	featureList=[]
# 	print "create polygons"
# 	itter=0
# 	for coordPair in polygons:
# 		for pt in coordPair:
# 			point.X=pt[0]
# 			point.Y=pt[1]
# 			array.add(point) 
# 		
# 		array.add(array.getObject(0))     
# 		polygon = arcpy.Polygon(array)
# 		array.removeAll()
# 		featureList.append(polygon)
# 	
# 		arcpy.CopyFeatures_management(featureList, tempPoly2)    
# 		arcpy.DefineProjection_management(tempPoly2, coordSys) 
# 		itter+=1
# 		print str(itter)+'/'+str(len(polygons))
# 	del point, array
	
	# initialize output part2
	meanOut=[]
	minOut=[]
	maxOut=[]
	print 'Calculate the mean, min, max and distance for each part of the swath'	
	ZSarea='tempTable.dbf'
	outZSaT=ZonalStatisticsAsTable(tempPoly2, 'FID', inRast, ZSarea)

	trows=arcpy.SearchCursor(ZSarea)
	for trow in trows:
		# save the 3 d area into a dictionary with the key as the FID, this is not in the loop
		# since these zonal statistics only need to be done once. 
		meanOut.append(trow.getValue('MEAN'))
		minOut.append(trow.getValue('MIN'))
		maxOut.append(trow.getValue('MAX'))
	del outZSaT, trows, trow        

	print 'delete the temporary files' 
	# arcpy.DeleteFeatures_management (tempPoly) 
	# arcpy.DeleteFeatures_management (tempPoly2) 
	# arcpy.DeleteFeatures_management (ZSAREA) 

	print 'Save Values'
	saveOut=[np.asarray(xOut), np.asarray(yOut), np.asarray(distOut),np.asarray(meanOut), np.asarray(minOut), np.asarray(maxOut)] 
	np.savetxt(fOutText, np.asarray(saveOut).T, delimiter=',')

	print "Make Plot"
	fig1=plt.figure(num=None, figsize=(3.35, 5), dpi=300, facecolor='w', edgecolor='w')
	#fig1.patch.set_alpha(0.0)

	ax = fig1.add_subplot(1,1,1)

	p1,=plt.plot(distOut, meanOut, lw=2, color='k', ls='-')
	p2,=plt.plot(distOut, minOut, lw=1, color='grey', ls='-')
	p3,=plt.plot(distOut, maxOut, lw=1, color='grey', ls='--')

	lg1=plt.legend([p1, p2, p3], ['Mean', 'Minimum', 'Maximum'], loc='upper left', fancybox=True,bbox_to_anchor=(0, 1))
	lg1.get_frame().set_edgecolor('grey')
	lg1.get_frame().set_linewidth(0.5)

	plt.tick_params(axis='both', which='major', labelsize=6)
	plt.tick_params(axis='both', which='minor', labelsize=6)

	plt.xlabel('Distance along Profile',fontsize=8)    
	plt.ylabel('Value',fontsize=8) 
	plt.title('Swath' ,fontsize=10)

	plt.savefig(fOutFigure, format='pdf')
	plt.show()
