print "importing modules"    
import arcpy
from arcpy import env
import matplotlib.pyplot as plt
import numpy as np
from krb_vecTools import *
env.overwriteOutput = 1
arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

####################### User Defined Inputs ##################################
env.workspace = r"C:\Users\Lab User\Billy\swath_profiler"
inPoly='xg1.shp'
#inRast='13jul15_dem.tif' # digital elevation map generated from july 2013 wv imagery
#inRast='julJul2013_corr3_dailyDisp_forArc.tif' # summer 2013 displacement (older version)
#inRast = 'summer_spring2013_difference.tif' # Difference map between summer and spring 2013 velocities
#inRast = '2013_dailyDisp_corrected_25feb14.tif' # june-july (summer) 2013 displacement (newer version)
#inRast = 'summer_spring2013_slidePct.tif' # Map of sliding as a percent of summer velocity
#inRast = 'marApr_2013_dailyDisp_11mar14.tif' # spring 2013 displacement
inRast = '01aug_27aug2013_totalDisp_corrected_forArc.tif' # 01 aug to 27 aug daily displacement
width=200.0 # this is in the units of the projection. 
height=10.0

fOutText='test.csv'
fOutFigure='test.pdf'

tempPoly='temp.shp' # interpolated polyline
tempPoly2='kenn_xg_bins.shp' # polygons for zonal statistics
ZSarea='tempTable.dbf' # zonal statistics table

####################### End User Defined Inputs ##############################

#### KRB August 2014 update notes
# 1. add something that checks the inputs (e.g. inPoly must be a line)





print "interpolate line to the raster"    
coordSys=arcpy.Describe(inPoly).spatialReference.exporttostring()
arcpy.InterpolateShape_3d (inRast, inPoly, tempPoly,height)

# initialize output part 1
xOut=[]
yOut=[]
print "get all points along the swath line"
ptArray=[]
rows=arcpy.SearchCursor(tempPoly)
for row in rows:
    points=row.shape.getPart(0)
    for point in points:
        xOut.append(point.X)
        yOut.append(point.Y)
        ptArray.append([point.X, point.Y])
del row, rows
     
print "calculate the slope and perpendicular slope along the line"
hwin=5
slope, perpSlope=calcSlopeAndPerpSlope(ptArray, hwin)

print "store the coordiates of the swath cross section lines"
polygons=[]
for i in range(len(ptArray)):
    polygons.append(makePoly(ptArray[i], perpSlope[i], width, height))

#calculate distance along the line
distOut=distanceAlongLine(xOut,yOut)

# initialize spatial containers
point=arcpy.Point()
array=arcpy.Array()
featureList=[]
print "create polygons"
itter=0
for coordList in polygons:
    for pt in coordList:
        point.X=pt[0]
        point.Y=pt[1]
        array.add(point) 
        
    array.add(array.getObject(0))     
    polygon = arcpy.Polygon(array)
    array.removeAll()
    featureList.append(polygon)
    
   
    itter+=1
    print str(itter)+'/'+str(len(polygons))
del point, array

arcpy.CopyFeatures_management(featureList, tempPoly2)    
arcpy.DefineProjection_management(tempPoly2, coordSys) 
	
# initialize output part2
meanOut=[]
minOut=[]
maxOut=[]
stdOut=[]
print 'Calculate the mean, min, max and distance for each part of the swath'	
outZSaT=ZonalStatisticsAsTable(tempPoly2, 'FID', inRast, ZSarea)

trows=arcpy.SearchCursor(ZSarea)
for trow in trows:
    # save the 3 d area into a dictionary with the key as the FID, this is not in the loop
    # since these zonal statistics only need to be done once. 
    meanOut.append(trow.getValue('MEAN'))
    minOut.append(trow.getValue('MIN'))
    maxOut.append(trow.getValue('MAX'))
    stdOut.append(trow.getValue('STD'))

del outZSaT, trows, trow        

print 'delete the temporary files' 
# arcpy.DeleteFeatures_management (tempPoly) 
# arcpy.DeleteFeatures_management (tempPoly2) 
# arcpy.DeleteFeatures_management (ZSAREA) 

print 'Save Values'
saveOut=[np.asarray(xOut), np.asarray(yOut), np.asarray(distOut),np.asarray(meanOut), np.asarray(minOut), np.asarray(maxOut)] 
np.savetxt(fOutText, np.asarray(saveOut).T, delimiter=',')


print 'Calculate the Location of Minimums and Maximums'
meanOut=np.array(meanOut)
minOut=np.array(minOut)
maxOut=np.array(maxOut)
distOut=np.array(distOut)
yOut=np.array(yOut)
xOut=np.array(xOut)


###Graph smoothed data###
plt.plot(meanOut,"k.-",label="original meanOut",alpha=.3)
plt.plot(smooth(meanOut,window_len=1, window='flat'),"-",label="smoothed d=1")
plt.plot(smooth(meanOut,window_len=3, window='flat'),"-",label="smoothed d=3")
plt.plot(smooth(meanOut,window_len=7, window='flat'),"-",label="smoothed d=5")
plt.title("Smoothing Window")
plt.grid(alpha=.3)
#plt.axis([0,2000,50,300])
plt.legend()
plt.show()

#differences=np.diff(smooth(meanOut,window_len=5, window='flat'))
#signChange=(differences[:-1])*(differences[1:])
#curvature=np.diff(differences)

## now try some different metrics for selecting ogives. Depending on what we are looking for we may want minima, maxima, etc. Can use slope, curvature, and sign change...
#
##selectedPoints=np.where(curvature>0)
#selectedPoints=np.where(signChange<0) #select indices for points where the sign changes, these are both the minima and the maximum
#selectCurvature=curvature[selectedPoints[0]]
#
#selInd=selectedPoints[0]
#
#
#
#curvIndTemp=np.where(selectCurvature>0)
#curvInd=curvIndTemp[0]
#
#finalIndices=selInd[curvInd]
#
#selectedPoint=finalIndices-2 #get indices, we must add one
#
#ind=np.array(range(curvature.size))
#plt.plot(ind, curvature, '.-')
#plt.plot(ind, signChange, '.-')
#plt.plot(ind[finalIndices], signChange[finalIndices], '.',color='r')
#plt.plot(ind[finalIndices], curvature[finalIndices], '.',color='r')
#plt.show()
#
#print 'Turn those back into Points and save to a Shapefile'
#plt.plot(xOut,yOut)
#plt.plot(xOut[selectedPoint],yOut[selectedPoint], '.', color='r')
#plt.show()
#
## A list of coordinate pairs
##
#selPointsOut='mySelectedPoints_r_p1.shp'
#
#xSelect=xOut[selectedPoint]
#ySelect=yOut[selectedPoint]
#distSelect=distOut[selectedPoint]
#
## Create an empty Point object
#point = arcpy.Point()
#
## A list to hold the PointGeometry objects
#pointGeometryList = []
#
## For each coordinate pair, populate the Point object and create
##  a new PointGeometry
#for i in range(ySelect.size):
#    point.X = xSelect[i]
#    point.Y = ySelect[i]
#    pointGeometry = arcpy.PointGeometry(point)
#    pointGeometryList.append(pointGeometry)
#
## Create a copy of the PointGeometry objects, by using pointGeometryList
##  as input to the CopyFeatures tool.
##
#arcpy.CopyFeatures_management(pointGeometryList, selPointsOut)   
#arcpy.DefineProjection_management(selPointsOut, coordSys) 
#
#arcpy.AddField_management (selPointsOut, 'MEAN', 'FLOAT')
#rows = arcpy.UpdateCursor(selPointsOut) 
#itter=0
#for row in rows:
#    # Fields from the table can be dynamically accessed from the row object.
#    #   Here fields named BUFFER_DISTANCE and ROAD_TYPE are used
#    row.MEAN = meanOut[itter]
#    rows.updateRow(row) 
#    itter+=1
#
## Delete cursor and row objects to remove locks on the data 
## 
#del row 
#del rows
#
#print "Calculate Ogive Legnths"
##Find the length of each ogive
#ogiveLength = np.diff(distSelect)
##Plot ogive lengths
#
#plt.plot(np.array(distSelect[1:73]), np.diff(distSelect), '.-', color='b')
#plt.title('Ogive Lengths', fontsize=10)
#plt.xlabel('Distance Along Profile', fontsize=8)
#plt.ylabel('Ogive Length', fontsize=8)
#plt.show()
#
#
#
#    
#
#print "Make Plot"
#fig1=plt.figure(num=None, figsize=(3.35, 5), dpi=300, facecolor='w', edgecolor='w')
##fig1.patch.set_alpha(0.0)
#
#ax = fig1.add_subplot(1,1,1)
#
#p1,=plt.plot(distOut, meanOut, lw=2, color='k', ls='-')
#plt.plot(np.array(distOut), np.array(meanOut), '.', color='b')
#plt.plot(np.array(distOut[selectedPoint]), np.array(meanOut[selectedPoint]), '.', color='r')
#
#p2,=plt.plot(distOut, minOut, lw=1, color='grey', ls='-')
#p3,=plt.plot(distOut, maxOut, lw=1, color='grey', ls='-')
#
#
#lg1=plt.legend([p1, p2, p3], ['Mean', 'Minimum', 'Maximum'], loc='upper left', fancybox=True,bbox_to_anchor=(0, 1))
#lg1.get_frame().set_edgecolor('grey')
#lg1.get_frame().set_linewidth(0.5)
#
#plt.tick_params(axis='both', which='major', labelsize=6)
#plt.tick_params(axis='both', which='minor', labelsize=6)
#
#plt.xlabel('Distance along Profile',fontsize=8)    
#plt.ylabel('Value',fontsize=8) 
#plt.title('Swath' ,fontsize=10)
#
#plt.savefig(fOutFigure, format='pdf')
plt.show()