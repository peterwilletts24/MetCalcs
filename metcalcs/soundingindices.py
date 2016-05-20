'''
Some sounding instability indicesbased on
http://weather.uwyo.edu/upperair/indices.html
Not tested really!
'''

from numpy import load, array, radians, sin, cos, linspace, mean, log, isnan, nan, nanmin, nanmax, nanmean, abs, zeros, exp, where,\
                  concatenate, diff

from metcalcs.soundingfuncs import *
import re

def KIndex(tempc, dew_tempc, pres):
    """
    """
    T850c = TempC850mb(tempc, pres)
    T500c = TempC500mb(tempc, pres)
    TD850c = DewTempC850mb(dew_tempc,pres)
    T700c = TempC700mb(tempc,pres)
    TD700c = DewTempC700mb(dew_tempc,pres)
    
    return (T850c - T500c) + TD850c - (T700c - TD700c)

def CrossTotalsIndex(tempc,dew_tempc, pres):
    """
    """

    TD850c = DewTempC850mb(dew_tempc,pres)
    T500c = TempC500mb(tempc, pres)
    
    return TD850c - T500c

def VerticalTotalsIndex(tempc,pres):
    """
    """
    T850c = TempC850mb(tempc,pres)
    T500c = TempC500mb(tempc, pres)
    
    return T850c - T500c

def TotalTotalsIndex(tempc, dew_tempc, pres):
    """
    """

    T850c = TempC850mb(tempc,pres)
    T500c = TempC500mb(tempc, pres)  
    TD850c = DewTempC850mb(tempc,pres)
    
    return (T850c - T500c) + (TD850c - T500c)

def LiftedIndex(tempc, dew_tempc, pres, heights, st_height):
    """
    LIFT	= T500 - Tparcel
		T500	= temperature in Celsius of the environment at 500 mb
		Tparcel	= 500 mb temperature in Celsius of a lifted parcel with 
        the average pressure, temperature, and dewpoint of the layer 500 m 
        above the surface 

    INPUTS: 
    500mb temp of parcel lifted from 850mb (C)
    500mb temp (C)
    pref: 

    OUTPUTS: Temp(C)

    Source: http://glossary.ametsoc.org/wiki/Stability_index
            http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.

    """
    T500c = TempC500mb(tempc, pres)
    
    
    t500c_lift_from_first_500m = TCParcelLiftedFromFirst500mTo500(tempc, dew_tempc, pres, heights, st_height)
    
    #print t500c_lift_from_first_500m
    #print T500c
        
    lift = T500c - t500c_lift_from_first_500m
    
    return lift

def ShowalterIndex(tempc, dew_tempc, pres):
    """
    SHOW	= T500 - Tparcel
		T500	= Temperature in Celsius at 500 mb
		Tparcel	= Temperature in Celsius at 500 mb of a parcel lifted from 850 mb 

    INPUTS: 
    500mb temp of parcel lifted from 850mb (C)
    500mb temp (C)
    pref: 

    OUTPUTS: Temp(C)

    Source: http://glossary.ametsoc.org/wiki/Stability_index
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.

    """
    t500c = TempC500mb(tempc,pres)
    
    t500c_lift_from_850 = TCParcelLiftedFrom850To500(tempc, dew_tempc, pres)
        
    show = t500c - t500c_lift_from_850
    
    return show
    
