# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 10:40:08 2012

@author: Katy
"""
import numpy as np
def pythagoreanTheorum(x,y):
    
    xOut=np.asarray(x)
    yOut=np.asarray(y)
    
    dist=(((xOut-xOut[0])**2)+(yOut-yOut[0])**2)**(0.5)
    return dist

def distanceAlongLine(x,y):   
    xOut=np.asarray(x)
    yOut=np.asarray(y)    
    segments=(((xOut[1:]-xOut[0:-1])**2)+(yOut[1:]-yOut[0:-1])**2)**(0.5)
    dist=np.cumsum(segments)
    dist=np.insert(dist, 0, 0)
    return dist
    
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
    
def makePoly(centerPt, slope, width, height, slopeUnits="radians"):
    if slopeUnits=="degrees":
        # convert to radians
        slope=slope/360*2*np.pi
    
    dx=np.cos(slope+(np.pi/2.0))*(height*0.5)
    dy=np.sin(slope+(np.pi/2.0))*(height*0.5)
    
    centerUp=[centerPt[0]+dx, centerPt[1]+dy]  
    centerDw=[centerPt[0]-dx, centerPt[1]-dy]
    
    dx2=np.cos(slope)*(width*0.5)
    dy2=np.sin(slope)*(width*0.5)
    
    NE=[centerUp[0]+dx2, centerUp[1]+dy2]
    NW=[centerUp[0]-dx2, centerUp[1]-dy2]
    
    SE=[centerDw[0]+dx2, centerDw[1]+dy2]
    SW=[centerDw[0]-dx2, centerDw[1]-dy2]    
    
    poly=[NW, NE, SE, SW]
    
    return poly 
    
def smoothTriangle(meanOut,degree,dropVals=False):        
    """performs moving triangle smoothing with a variable degree."""        
    """note that if dropVals is False, output length will be identical        
    to input length, but with copies of data at the flanking regions"""        
    triangle=np.array(range(degree)+[degree]+range(degree)[::-1])+1        
    smoothed=[]        
    for i in range(degree,len(meanOut)-degree*2):                
        point=meanOut[i:i+len(triangle)]*triangle                
        smoothed.append(sum(point)/sum(triangle))        
    if dropVals: return smoothed        
    smoothed=[smoothed[0]]*(degree+degree/2)+smoothed        
    while len(smoothed)<len(meanOut):smoothed.append(smoothed[-1])        
    return smoothed

def smoothKaiser(x,beta):
    """ kaiser window smoothing """
    window_len=11
    # extending the data at beginning and at the end
    # to apply the window at the borders
    s = np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    w = np.kaiser(window_len,beta)
    y = np.convolve(w/w.sum(),s,mode='valid')
    return y[5:len(y)-5]
 
 
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y
    
def momCoolFxn(randomValue):
        #if you wanted to do something with randomValue you'd do it here
        return 'myMomIsAwesome'