print "importing modules"    
import numpy as np
from vecTools import *

from osgeo import ogr # modify code to just use ogr, not arcGIS

env.workspace = r"C:\Users\Lab User\Desktop\katy\swathProfile\data"
inPoly='crossSectionLine.shp'
inRast='orthoWV01_09JUL132056451-P1BS-1020010008B20800_u08mm26907.tif'
width=100.0 # this is in the units of the projection. 
height=5.0

fOutText='fOutCSV.csv'
fOutFigure='fOut.pdf'

tempPoly='temp.shp'
tempPoly2='temp2.shp'
tempPoly3='temp3.shp'

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
for coordPair in polygons:
    for pt in coordPair:
        point.X=pt[0]
        point.Y=pt[1]
        array.add(point) 
        
    array.add(array.getObject(0))     
    polygon = arcpy.Polygon(array)
    array.removeAll()
    featureList.append(polygon)
    
    arcpy.CopyFeatures_management(featureList, tempPoly2)    
    arcpy.DefineProjection_management(tempPoly2, coordSys) 
    itter+=1
    print str(itter)+'/'+str(len(polygons))
del point, array
	
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
