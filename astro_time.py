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
    
def ymd2doy(year, month, day):
    mos = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (is_leap_year(year)):
        mos[1] = 29

    idx = month - 1
    
    doy = np.sum(mos[:idx]) + day
    return doy
    
def doy2ymd(day_of_year, year):
    mos = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (is_leap_year(year)):
        mos[1] = 29
    
    temp = 0
    idx = 0
    while (temp < day_of_year):
        temp += mos[idx]
        
        if(temp >= day_of_year):
            month = idx + 1
            if (month == 1):
                day = day_of_year
            else:
                day = day_of_year - np.sum(mos[:idx])
            return (month, day)
        
        idx += 1
        
def ymdhms2days(year, month, day, hour, minutes, seconds):
    doy = ymd2doy(year, month, day)
    days = doy + hour/24.0 + minutes/1440.0 + seconds/86400.0
    return days
    
def days2ymdhms(days, year):
    doy = np.trunc(days)
    (month, day) = doy2ymd(doy, year)
    
    tau = (days - doy)*86400.0
    (hours, minutes, seconds) = time2hms(tau)
    
    return (month, day, hours, minutes, seconds)
    
def jd2gregorian(jd):
    t1900 = (jd - 2415019.5)/365.25
    year = 1900 + np.trunc(t1900)
    leap_years = np.trunc((year - 1900 - 1)*0.25)
    days = (jd - 2415019.5) - ((year - 1900)*365.0 + leap_years)
    
    if days < 1.0:
        year = year - 1
        leap_years = np.trunc((year - 1900 - 1)*0.25)
        days = (jd - 2415019.5) - ((year - 1900)*365.0 + leap_years)
        
    (month, day, hours, minutes, seconds) = days2ymdhms(days, year)
    return (year, month, day, hours, minutes, seconds)