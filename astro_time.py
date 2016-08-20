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