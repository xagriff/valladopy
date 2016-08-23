# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 17:34:45 2016

@author: Alex
"""

import numpy as np
import math

def julian_date(yr, mo, d, hr, minute, sec, leap_sec=False):
    x = (7*(yr + np.trunc((mo + 9)/12)))/4.0
    y = (275*mo)/9.0
    if leap_sec:
        t = 61.0
    else:
        t = 60.0
    z = (sec/t + minute)/60.0 + hr
    jd = 367*yr - np.trunc(x) + np.trunc(y) + d + 1721013.5 + z/24.0
    return jd
    
def find_gmst(jd_ut1):
    t_ut1 = (jd_ut1 - 2451545.0)/36525.0
    theta_gmst = 67310.54841 + (876600*3600.0 + 8640184.812866)*t_ut1 + 0.093104*t_ut1*t_ut1 - 6.2E-6*t_ut1*t_ut1*t_ut1
    theta_gmst = math.fmod(theta_gmst, 86400.0)
    theta_gmst /= 240.0
    if theta_gmst < 0:
        theta_gmst += 360.0
    return theta_gmst
    
def find_lst(theta_gmst, lon):
    theta_lst = theta_gmst + lon
    return theta_lst
    
def dms2rad(degrees, minutes, seconds):
    rad = (degrees + minutes/60.0 + seconds/3600.0)*(math.pi/180.0)
    return rad
    
def rad2dms(rad):
    temp = rad*(180.0/math.pi)
    degrees = np.trunc(temp)
    minutes = np.trunc((temp - degrees)*60.0)
    seconds = (temp - degrees - minutes/60.0)*3600.0
    return (degrees, minutes, seconds)
    
def hms2rad(hours, minutes, seconds):
    rad = 15*(hours + minutes/60.0 + seconds/3600.0)*math.pi/180.0
    return rad
    
def rad2hms(rad):
    temp = rad*180.0/(15.0*math.pi)
    hours = np.trunc(temp)
    minutes = np.trunc((temp - hours)*60.0)
    seconds = (temp - hours - minutes/60.0)*3600.0
    return (hours, minutes, seconds)
    
def is_leap_year(yr):
    if (np.remainder(yr,4) != 0):
        return False
    else:
        if (np.remainder(yr, 100) == 0):
            if (np.remainder(yr, 400) == 0):
                return True
            else:
                return False
        else:
            return True
            
def time2hms(time_in_seconds):
    temp = time_in_seconds/3600.0
    hours = np.trunc(temp)
    minutes = np.trunc((temp - hours)*60)
    seconds = (temp - hours - minutes/60)*3600
    
    return (hours, minutes, seconds)
    
def hms2time(hours, minutes, seconds):
    tau = 3600.0*hours + 60.0*minutes + seconds
    return tau