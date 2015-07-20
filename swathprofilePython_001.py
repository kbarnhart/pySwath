print "importing modules"    


# ultimately we want to make this a function, with the following input information:

### USER DEFINED VARIABLES ####

folderPath = u'/Users/wiar9509/git/pySwath'
inPoly='/test/testMultiLine.shp' # input shapefile filename along which raster will be sampled
inRast='vv.tif' # input raster filename to sample
width=100.0 # across-flow box dimension [this is in the units of the projection]
height=5.0 # along-flow box dimension

fOutText='fOutCSV.csv' # filename for textfile output
fOutFigure='fOut.pdf' # filename for figure output

### LOADING MODULES ###

import os
os.chdir(folderPath)
import numpy as np
from krb_vecTools import *
import matplotlib.pyplot as plt

from osgeo import ogr # modify code to just use ogr, not arcGIS


### PROCESSING ###

tempPoly='temp.shp'
tempPoly2='temp2.shp'
tempPoly3='temp3.shp'

print "Opening shapefile with OGR"
shapefile=folderPath+inPoly
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(shapefile, 0)
layer = dataSource.GetLayer()

# get all points along shapefile
numLines=layer.GetFeatureCount()

for lineNo in range(numLines): # need to add additional info for output files if there are more than one line. 20jul15 WHA - did this in 'name' variable below
	
	lineX=[]
	lineY=[]
 	
 	feat=layer.GetFeature(lineNo) # Highlights current line
	geom=feat.geometry()
	name=feat.GetField(1) # Gets current line name
	
	print "Getting coordinates along line: "+name
	
	numPoints=geom.GetPointCount()
	
	for nP in range(numPoints): # Iterates over points and appends x,y (easting, northing) coords from current line
		lineX.append(geom.GetPoint(nP)[0])
		lineY.append(geom.GetPoint(nP)[1])
		
	
		
	# interpolate to the spacing specified by height (get x,y coordinates along the line at a distance of height apart)
	# put them into an array called ptArray, that is a list of lists, with each interior list being [pointX, pointY]
	# e.g. ptArray.append([pointX, pointY]) 
	
	
   
# 	coordSys=arcpy.Describe(inPoly).spatialReference.exporttostring()
# 	arcpy.InterpolateShape_3d (inRast, inPoly, tempPoly,height)
# 
# 	
	# xOut=[]
# 	yOut=[]
# 	print "get all points along the swath line"
# 	ptArray=[]
# 
# 
# 
# 
# 	rows=arcpy.SearchCursor(tempPoly)
# 	for row in rows:
# 		points=row.shape.getPart(0)
# 		for point in points:
# 			xOut.append(point.X)
# 			yOut.append(point.Y)
# 			ptArray.append([point.X, point.Y])
# 	del row, rows
	 
	print "calculate the slope and perpendicular slope along the line"
	hwin=5
	slope, perpSlope=calcSlopeAndPerpSlope(ptArray, hwin)

	print "store the coordiates of the swath cross section lines"
	polygons=[]
	for i in range(len(ptArray)):
		polygons.append(makePoly(ptArray[i], perpSlope[i], width, height))

	#calculate distance along the line
	distOut=distanceAlongLine(xOut,yOut)

	# initialize spatial containers (Need to revise with OGR output)
	
	
	
	
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
