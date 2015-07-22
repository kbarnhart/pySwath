# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 10:40:08 2012

@author: Katy
"""
import numpy as np

def calculateSlope(pt1,pt2, units="radians"):
    
    ''' calculates slope between horizontal and line defined by p1=[x1,y1] and 
    p2=[x2,y2] at the point defined by [x1, y1] using np.arctan2'''
    
    if (np.asarray(pt1).shape[0]==np.asarray(pt2).shape[0]) and (np.asarray(pt1).size%2==0 ) and(np.asarray(pt2).size%2==0 ):   
    
        if np.asarray(pt1).size>2: 
            # reqires unzipping...
            x1, y1=zip(*pt1)    
            x2, y2=zip(*pt2)
            x=[]
            y=[]
            # f-you tuples...
            for i in range(len(x1)):
                x.append(x2[i]-x1[i])
                y.append(y2[i]-y1[i])                       
        else:
            x=pt2[0]-pt1[0]
            y=pt2[1]-pt1[1]  
            
        slope=np.arctan2(y,x)    
        if units=="degrees":
            # convert to degrees
            slope=slope/(2*np.pi)*360        
        return slope

def calcSlopeAndPerpSlope(ptArray, hwin, units="radians"):

    slope=[]
    perpSlope=[]
    
    for i in range(hwin): 
        firstPt=ptArray[0]
        lastPt=ptArray[i+hwin]
        slope.append(calculateSlope(firstPt, lastPt,units))
        perpSlope.append(((slope[i]+2*np.pi)-(np.pi/2))%(2*np.pi))
        
    for i in range(hwin,len(ptArray)-hwin):
        firstPt=ptArray[i-hwin]
        lastPt=ptArray[i+hwin]
        slope.append(calculateSlope(firstPt, lastPt,units))
        perpSlope.append(((slope[i]+2*np.pi)-(np.pi/2))%(2*np.pi))
        
    for i in range(len(ptArray)-hwin, len(ptArray)): 
        firstPt=ptArray[i-hwin]
        lastPt=ptArray[-1]
        slope.append(calculateSlope(firstPt, lastPt,units))
        perpSlope.append(((slope[i]+2*np.pi)-(np.pi/2))%(2*np.pi))    
    return slope, perpSlope

def makeLine(centerPt, slope, halfLength, slopeUnits="radians"):
    if slopeUnits=="degrees":
        # convert to radians
        slope=slope/360*2*np.pi
    
    dx=np.cos(slope)*halfLength
    dy=np.sin(slope)*halfLength    
    firstPt=[centerPt[0]-dx, centerPt[1]-dy]
    lastPt=[centerPt[0]+dx, centerPt[1]+dy]
    line=[firstPt, centerPt, lastPt]    
    return line