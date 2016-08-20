# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 17:34:45 2016

@author: Alex
"""

import numpy as np

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