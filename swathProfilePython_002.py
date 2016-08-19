### Swath profiler code ###
# Written by William Armstrong and Katy Barnhart, with code modified from KRB's Arc-based profiler
# Begun 15 July 2015
# Finished (sort of) 23 July 2015

###### LOADING MODULES ######

import os
folderPath='/Users/wiar9509/Documents/generalScripts/swath/' # where are swath profiler and vecTools located?
os.chdir(folderPath)
import sys
import numpy as np
from krb_vecTools import *
import scipy as sci 
import grass.script as grass
from osgeo import ogr
from osgeo import osr
from osgeo import gdal
from gdalconst import *
import csv
import glob

def swathProfiler(imagePath,inRast,inPoly):

	print "STARTING SWATH PROFILER"
	print "Importing modules"    


	##### USER DEFINED VARIABLES ######

	# ultimately we want to make this a function, with the following input information:
	#folderPath='/Users/wiar9509/git/pySwath/'
	os.chdir(folderPath)
	#inPoly='/test/testMultiLine.shp' # input shapefile filename along which raster will be sampled
	#inPoly='test9.shp' # input shapefile filename along which raster will be sampled
	#inPoly='wrangellsTransectsUtm7n.shp' # input shapefile filename along which raster will be sampled
	#='stEliasTransectsUtm7n.shp' # input shapefile filename along which raster will be sampled
	#inRast='vv.tif' # input raster filename to sample
	#inRast='LC80640172013166LGN0_LC80640172014137LGN0_2013-334_336_20_69_hp_filt_3.0_vv.tif'
	#inRast='LC80640172013166LGN0_LC80640172013214LGN0_2013-190_48_20_69_hp_filt_3.0_vv.tif'
	#inRast='wrangells_gdem.tif'
	width=750.0 # across-flow box dimension [this is in the units of the projection]
	height=500.0 # along-flow box dimension

	###### PROCESSING ######
	print 'PROFILING RASTER: '+inRast
	os.chdir(folderPath)
	# Prepare outputs
	outTextSuffix='stats.csv' # filename suffix for textfile output
	figOutSuffix='statsFig.pdf' # filename suffix for figure output
	#inRastOut=inRast[0:len(inRast)-4] # trim the .tif from filename
	rastSplit=inRast.split('_') # break up input raster string
	masterImage=rastSplit[0]
	slaveImage=rastSplit[1]

	# Create output folder structure if it doesn't exist
	if not os.path.exists(imagePath+'swathOutputs/'):
		os.mkdir(imagePath+'swathOutputs')
		os.mkdir(imagePath+'swathOutputs/textfiles')
		os.mkdir(imagePath+'swathOutputs/figures')
		os.mkdir(imagePath+'swathOutputs/shapefiles')

	# Create the spatial reference
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(32607) # This hard-coded to UTM 7N w/ WGS84 datum. Would be nice to have this defined based of input shapefile

	print "Opening shapefile with OGR"
	shapefile=folderPath+'shapefileInputs/'+inPoly
	driver = ogr.GetDriverByName("ESRI Shapefile")
	dataSource = driver.Open(shapefile, 1)
	layer = dataSource.GetLayer()
	crs = layer.GetSpatialRef() # Coordinate reference system
	
	# Function to find raster bounding box from : https://www.siafoo.net/snippet/69/rev/1
	def findGDALCoordinates(path):
		if not os.path.isfile(path):
			return []
		data = gdal.Open(path,GA_ReadOnly)
		if data is None:
			return []
		geoTransform = data.GetGeoTransform()
		minx = geoTransform[0]
		maxy = geoTransform[3]
		maxx = minx + geoTransform[1]*data.RasterXSize
		miny = maxy + geoTransform[5]*data.RasterYSize
		return [minx,miny,maxx,maxy]
	
	# Corner coordinates for input raster. Was going to add a test to not sample if vector not within extent of raster, but the no-data areas complicate things.
	[rMinX,rMinY,rMaxX,rMaxY] = findGDALCoordinates(imagePath+inRast)
	
	# Changing filepath so GRASS won't balk. OGR needs it a different way it seems. Just need to trim the starting 'u' before filepath
	if not folderPath[0]=='/' and folderPath[1]=='/':
		pathLen=len(folderPath)
		grassPath=folderPath[1:pathLen]
	else:
		grassPath=folderPath

	print "Checking projection compatibility between polylines and raster"
	# Read raster projection
	#rastProj=grass.parse_command("g.proj",flags='j',georef=grassPath+inRast)
	rastProj=grass.parse_command("g.proj",flags='j',georef=imagePath+inRast)
	# Read polyline projection
	#polyProj=grass.parse_command("g.proj",flags='j',georef=grassPath+inPoly)
	polyProj=grass.parse_command("g.proj",flags='j',georef=folderPath+'shapefileInputs/'+inPoly)

	if rastProj.values()[4]==polyProj.values()[4] and rastProj.values()[6]==polyProj.values()[6]:
		print 'Raster and polylines have matching projection'
# 	elif rastProj.values()[2]=='longlat': # If raster is wgs84 will convert to utm 7n/wgs84 (abandoned this on 24Sep2015 16:02; not currently working).
# 		# Input projection is WGS84
# 		inSpatialRef = osr.SpatialReference()
# 		inSpatialRef.ImportFromEPSG(4326)
# 		# Output projection is UTM 7N/WGS 84
# 		outSpatialRef = osr.SpatialReference()
# 		outSpatialRef.ImportFromEPSG(32607) # Not sure if this number right, copied for elsewhere
# 		coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef) # transformation matrix
	else:
		print 'Raster and polylines DO NOT have matching projection, terminating script'
		#sys.exit(1)

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
		extent=geom.GetEnvelope() # Bounding box of profile vector (could be used for testing)
		
		# Get profile extent
		pMinX=extent[0] # profile min x
		pMaxX=extent[1] # profile max x
		pMinY=extent[2] # profile min y
		pMaxY=extent[3] # profile max y
		
		profileExtent=[]
		profileExtent.append([pMinX,pMinY,pMaxX,pMaxY])	
		
		rasterExtent=[]
		rasterExtent.append([rMinX,rMinY,rMaxX,rMaxY])
				
		# Test if profile polyline is contained within raster. Right now needs to be entirely in shapefile. Need to modify later.
		def lineWithinPolygon(profileExtent,rasterExtent):
			# If line completely contained within polygon
			if profileExtent[0][0] >= rasterExtent[0][0] and profileExtent[0][2] <= rasterExtent[0][2] and profileExtent[0][1] >= rasterExtent[0][1] and profileExtent[0][3] <= rasterExtent[0][3]  :
				value = 1
				return value
		
		# This '1' if profile line is fully contained within raster extent
		lineContained = lineWithinPolygon(profileExtent,rasterExtent)
		
		# Naming output
		if inRast[0:2] == 'LC':
			outTextFn=imagePath+'swathOutputs/textfiles/'+name+'_'+masterImage+'_'+slaveImage+'_'+outTextSuffix
			rastName=masterImage+'_'+slaveImage
		else:
			rastName=inRast
	
		rastNameLen=len(inRast)

		if rastName[rastNameLen-4:rastNameLen] == '.tif':
			rastName=rastName[0:rastNameLen-4]
		outTextFn=imagePath+'swathOutputs/textfiles/'+name+'_'+rastName+'_'+outTextSuffix
	
		print "PROCESSING TRANSECT "+name.upper()
		print "TRANSECT "+str(lineNo+1)+" OF "+str(numLines)
		print "Getting coordinates along line"

		# Only run rest of script if transect line is within raster extent
		if lineContained == 1:
			print "Transect within raster extent"
			
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

			# Vectools to generate polygons (rectangles) for zonal stats 
			print "Calculating the slope and perpendicular slope along the line"
			hwin=5
			slope, perpSlope=calcSlopeAndPerpSlope(ptArray, hwin) # don't get this hwin thing

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
			if os.path.exists(folderPath+'/swathOutputs/shapefiles/'+name+"polygons.shp"):
				driver.DeleteDataSource(folderPath+'/swathOutputs/shapefiles/'+name+"polygons.shp") # error if data source already exists
			newDataSource=driver.CreateDataSource(folderPath+'/swathOutputs/shapefiles/'+name+"polygons.shp")
			newLayer=newDataSource.CreateLayer(name+"polygons",srs,geom_type=ogr.wkbPolygon)
			fieldDefn=ogr.FieldDefn('id',ogr.OFTInteger)
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
			grass.run_command("r.in.gdal",flags='e',overwrite=True,quiet=True,input=imagePath+inRast,output='sampleRast')
			# Read in polygons to use as sampling bins
			grass.run_command("v.in.ogr",overwrite=True,quiet=True,input=folderPath+'/swathOutputs/shapefiles/'+name+"polygons.shp",output='samplePolys')
			# Make sure computational region contains polygons
			grass.run_command("g.region",quiet=True,vec='samplePolys',res='10')
			# Convert polygons vector into 'zone' raster for sampling
			grass.run_command("v.to.rast",overwrite=True,quiet=True,input='samplePolys',output='zoneRast',use='attr',attr='id')
			# Run univariate statistics on sample raster, using zone raster to bin and populate new rows
			grass.run_command("r.univar",overwrite=True,quiet=True,flags='t',map='sampleRast', zones='zoneRast', out=outTextFn, sep=',',)
			# You now have a file named $outTextFilename in the folder in which this swath profiler code lives

			### END GRASS PORTION ###

			# Initializing stats
			min=[]
			max=[]
			rangeVals=[]
			mean=[]
			stdDev=[]
			sum=[]

			# Read in stats file
	# 		print 'Reading in statistics file and updating shapefile'
	# 		statsFile=np.genfromtxt(outTextFn,delimiter=',',skip_header=1)
	# 		statsLen=len(statsFile)

			# Below appends the along-transect coordinate to the stats textfile
			# http://stackoverflow.com/questions/23682236/add-a-new-column-to-an-existing-csv-file
			# This isn't working as of 27 July 2015
		# 	csvfile = outTextFn
		# 	fnLen=len(csvfile)
		# 	newName = outTextFn[0:fnLen-4]+'_new'+outTextFn[fnLen-4:fnLen]
		# 	with open(csvfile, 'r') as fin, open(newName, 'w') as fout:
		# 		reader = csv.reader(fin, delimiter=',', lineterminator='\n')
		# 		writer = csv.writer(fout, delimiter=',', lineterminator='\n')
		# 		writer.writerow(next(reader) + ['dist']
		# 		for row, val in zip(reader, lineDist)
		# 			writer.writerow(row + [lineDist])


			# Outputting distance along line textfile (this because can't get the add column above working)
			np.savetxt(imagePath+'swathOutputs/textfiles/'+name+'_'+rastName+'_alongTransectDist.csv',lineDist,delimiter=',')
			
		else:
			print "Transect not within raster extent; not processing"

	# #Opening layer to populate stats values in shapefile created above
	# 		updateDatasource = driver.Open(folderPath+'/swathOutputs/shapefiles/'+name+"polygons.shp", update=1)
	# 		updateLayer = updateDatasource.GetLayer()
	# 
	# 		for line in range(statsLen): # First value of these is nan
	# 
	# 			min.append(statsFile[line][4]) # minimum value of raster within polygon
	# 			max.append(statsFile[line][5]) # maximum value
	# 			rangeVals.append(statsFile[line][6]) # range of values
	# 			mean.append(statsFile[line][7]) # mean value
	# 			stdDev.append(statsFile[line][9]) # standard deviation of values
	# 			sum.append(statsFile[line][12])	# sum of values
	# 
	# 			# Iterating over boxes to set stats
	# 			feat=updateLayer.GetFeature(line)
	# 			feat.SetField('min',min[line])
	# 			feat.SetField('max',max[line])
	# 			feat.SetField('mean',mean[line])
	# 			feat.SetField('range',rangeVals[line])
	# 			feat.SetField('stddev',stdDev[line])
	# 			feat.SetField('sum',sum[line])
	# 
	# 			# Update the feature in the layer
	# 			updateLayer.SetFeature(feat)
	# 
	# 			# Destroy the feature to free resources
	# 			feat.Destroy()
	# 
	# #Destroy the data source to free resources
	# 		updateDatasource.Destroy()
	# 	
			# Loop to next line in shapefile
		

#### FOR LOOPING THIS SCRIPT THROUGH MULTIPLE IMAGES

imageryPath='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/corr/pycorr/vv_files/filtered'
path='064'
row='017'
#inPoly='wrangellsTransectsUtm7n.shp' # input shapefile filename along which raster will be sampled
#inPoly='stEliasTransectsUtm7n.shp' # input shapefile filename along which raster will be sampled
#inPoly='wStEtransectsUtm7n.shp' # input shapefile filename along which raster will be sampled
inPoly='newWStETransectsUtm7n_24sep15.shp'
#inRast='LC80630172013175LGN0_LC80630172014066LGN0_2013-303_256_24_40_10_hp_filt_3.0_vv.tif'
#masterDoy='2013166' # make [] to not specify
#imagePath='/Users/wiar9509/Documents/CU2014-2015/wrangellStElias/dtm/'
#inRast='wrangellStElias_gdem.tif'
#inRast='stelias_gdem.tif'

masterDoy=[]
#### PROCESSING

os.chdir(imageryPath)

if not masterDoy:
	fileList=glob.glob('LC8'+path+row+'*')
else:
	fileList=glob.glob('LC8'+path+row+masterDoy+'*')

fileNum=len(fileList)
iter = 0

for file in range(fileNum):
	fileNow=fileList[file]
	nameLen = len(fileNow)
	fileEnd = fileNow[nameLen-4:nameLen]
	if fileEnd == '.tif':
		swathProfiler(imageryPath,fileNow,inPoly)

	