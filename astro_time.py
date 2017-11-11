# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 17:34:45 2016

@author: Alex
"""

import numpy as np
import math

def julian_date(yr, mo, d, hr, minute, sec, leap_sec=False):
    """Converts a date & time to Julian Date
    
    Converts a date & time (yr, month, day, hour, second) to a Julian Date.
    Valid for any time system (UT1, UTC, AT, etc.) but should be identified to
    avoid confusion. For reference, see Algorithm 14 in Vallado (Fourth 
    Edition), Section 3.5 pg 185.
    
    Parameters
    ----------
    yr: int
        Four digit year
    mo: int
        Month
    d:  int
        Day (of month)
    hr: int
        Hour (24-hr based)
    minute: int
        Minute
    sec: double
        Seconds
    leap_sec: boolean, optional, default = False
        Flag if time is during leap second
    
    Returns
    -------
    jd: double
        Date/Time as a Julian Date
    """
    
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
    """Finds Greenwich Mean Sidereal Time given UT1
    
    Finds the Greenwich Mean Sidereal Time (GMST) for a supplied UT1 Julian 
    Date. For reference, see Algorithm 15 in Vallado (Fourth Edition), 
    Section 3.5 pg 188.
    
    Parameters
    ----------
    jd_ut1: double
        The UT1 Julian Date
        
    Returns
    -------
    theta_gmst: double
        The Greenwich Mean Sidereal Time, expressed as angle in degrees       
    """
    t_ut1 = (jd_ut1 - 2451545.0)/36525.0
    theta_gmst = 67310.54841 + (876600*3600.0 + 8640184.812866)*t_ut1 + 0.093104*t_ut1*t_ut1 - 6.2E-6*t_ut1*t_ut1*t_ut1
    theta_gmst = math.fmod(theta_gmst, 86400.0)
    theta_gmst /= 240.0
    if theta_gmst < 0:
        theta_gmst += 360.0
    return theta_gmst
    
def find_lst(theta_gmst, lon):
    """Finds the Local Sidereal Time given GMST and longitude
    
    Finds the Local Sidereal Time (LST) for a supplied GMST and Longitude. For
    reference, see Algorithm 15 in Vallado (Fourth Edition), Section 3.5 pg 188
    
    Parameters
    ----------
    theta_gmst: double
        GMST as an angle in degrees
    longitude: double
        Longitude of the site of interest, in degrees
        
    Returns
    -------
    theta_lst: double
        The Local Sidereal Time, expressed as an angle in degrees
    """
    
    theta_lst = theta_gmst + lon
    return theta_lst
    
def dms2rad(degrees, minutes, seconds):
    """ Converts degrees, minutes, seconds to radians.
    
    Converts degrees, minutes, seconds to radians. For reference, see Algorithm
    17 in Vallado (Fourth Edition), Section 3.5 pg 197
    
    Parameters
    ----------
    degrees: double
        degrees part of angle
    minutes: double
        minutes part of angle
    seconds: double
        seconds part of angle
        
    Returns
    -------
    rad: double
        angle in radians
    """
    rad = (degrees + minutes/60.0 + seconds/3600.0)*(math.pi/180.0)
    return rad
    
def rad2dms(rad):
    """Converts an angle in radians to Degrees, Minutes, Seconds
    
    Converts and angle in radians to Degrees, Minutes, Seconds. For reference, 
    see Algorithm 18 in Vallado (Fourth Edition), Section 3.5 pg 197
    
    Parameters
    ----------
    rad: double
        Angle in radians
        
    Returns
    -------
    (degrees, minutes, seconds): tuple
        Angle parts in degrees, minutes, and seconds
    """
    temp = rad*(180.0/math.pi)
    degrees = np.trunc(temp)
    minutes = np.trunc((temp - degrees)*60.0)
    seconds = (temp - degrees - minutes/60.0)*3600.0
    return (degrees, minutes, seconds)
    
def hms2rad(hours, minutes, seconds):
    """Converts a time (hours, minutes, and seconds) to an angle (radians).
    
    Converts a time (hours, minutes, and seconds) to an angle in radians. For
    reference, see Algorithm 19 in Vallado (Fourth Edition), Section 3.5 pg 198
    
    Parameters
    ----------
    hours: double
        hours portion of time
    minutes: double
        minutes portion of time
    seconds: double
        seconds portion of time
        
    Returns
    -------
    rad: double
        angle representation in radians
    """
    rad = 15*(hours + minutes/60.0 + seconds/3600.0)*math.pi/180.0
    return rad
    
def rad2hms(rad):
    """Converts an angle (in radians) to a time (in hours, minutes, seconds)
    
    Converts an angle (radians) to a time (hours, minutes, seconds). For
    refrence, see Algorithm 20 in Vallado (Fourth Edition), Section 3.5 pg 198
    
    Parameters
    ----------
    rad: double
        angle in radians
        
    Returns
    -------
    (hours, minutes, seconds): tuple
        time in hours, minutes, and seconds
    """
    temp = rad*180.0/(15.0*math.pi)
    hours = np.trunc(temp)
    minutes = np.trunc((temp - hours)*60.0)
    seconds = (temp - hours - minutes/60.0)*3600.0
    return (hours, minutes, seconds)
    
def is_leap_year(yr):
    """Determines if the year is a leap year
    
    Determines if the year is a leap year. For reference, see Section 3.6.4 of
    Vallado (Fourth Edition), pg 200
    
    Parameters
    ----------
    yr: int
        Four digit year
        
    Returns
    -------
    is_leap: boolean
        Flag indicating if the year is a leap year
    """
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
    """Computes the time in Hours, Minutes, Seconds from seconds into a day
    
    Computes the time in Hours, Minutes, and Seconds from the time expressed 
    as seconds into a day.  For reference, see Algorithm 21 in Vallado 
    (Fourth Edition), Section 3.6.3 pg 199
    
    Parameters
    ----------
    time_in_seconds: double
        time expessed as seconds into a day
        
    Returns
    -------
    (hours, minutes, seconds): tuple
        time expressed as hours, minutes, seconds
    """
    temp = time_in_seconds/3600.0
    hours = np.trunc(temp)
    minutes = np.trunc((temp - hours)*60)
    seconds = (temp - hours - minutes/60)*3600
    
    return (hours, minutes, seconds)
    
def hms2time(hours, minutes, seconds):
    """Computes the time in seconds into the day from hours, minutes, seconds
    
    Computes the time in seconds into the day from the time expressed as 
    hours, minutes, seconds. For reference, see Section 3.6.3 in Vallado 
    (Foiurth Edition), pg 199
    
    Parameters
    ----------
    hours: double
        hours part of the time
    minutes: double
        minutes part of the time
    seconds: double
        seconds part of the time
        
    Returns
    -------
    tau: double
        the time as seconds into the day
    """
    tau = 3600.0*hours + 60.0*minutes + seconds
    return tau
    
def ymd2doy(year, month, day):
    """Computes the day of the year from the year, month, and day
    
    Computes the day of the year from the year, month, and day of a date. For
    reference, see Section 3.6.4 of Vallado (Fourth Edition) pg 200
    
    Parameters
    ----------
    year: int
        year (needed for leap year test)
    month: int
        month
    day: int
        day
        
    Returns
    -------
    doy: int
        day of the year
    """       
    mos = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (is_leap_year(year)):
        mos[1] = 29

    idx = month - 1
    
    doy = np.sum(mos[:idx]) + day
    return doy
    
def doy2ymd(day_of_year, year):
    """Computes the month and day, given the year and day of the year
    
    Computes the month and day, given the year and day of the year. For 
    reference, see Section 3.6.4 in Vallado (Fourth Edition), pg 200
    
    Parameters
    ----------
    day_of_year: int
        the day of the year
    year: int
        the year (needed for leap year test)
        
    Returns
    -------
    (month, day): tuple
        the month and day of the date
    """
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
    """Computes the decimal day for a give date (y, m, d) and time (h, m, s)
    
    Computes the decimal day for a given date (years, months, days) and time
    (hours, minutes, seconds). For reference, see Section 3.6.5 in Vallado 
    (Fourth Edition) pg 201
    
    Parameters
    ----------
    year: int
        year (YYYY)
    month: int
        month
    day: int
        day
    hour: int
        hour
    minutes: int
        minutes
    seconds: double
        seconds
        
    Returns
    -------
    days: double
        decimal day (and partial day) of the year
    """
    doy = ymd2doy(year, month, day)
    days = doy + hour/24.0 + minutes/1440.0 + seconds/86400.0
    return days
    
def days2ymdhms(days, year):
    """Computes the month, day, hrs, mins, and secs from the year and decimal day
    
    Computes the month, day, hours, minutes, and seconds from the year and 
    decimal day of the year (including partial day). For reference, see 
    Section 3.6.5 in Vallado (Fourth Edition), pg 201
    
    Parameters
    ----------
    days: double
        decimal day (and partial day) of the year
    year: int
        year (YYYY)
        
    Returns
    -------
    (month, day, hours, minutes, seconds): tuple
        the month, day, hours, minutes, and seconds of the date/time
    """
    doy = np.trunc(days)
    (month, day) = doy2ymd(doy, year)
    
    tau = (days - doy)*86400.0
    (hours, minutes, seconds) = time2hms(tau)
    
    return (month, day, hours, minutes, seconds)
    
def jd2gregorian(jd):
    """Computes the components of the Gregorian Date (y, mo, d, h, m, s) from 
    a Julian Date
    
    Computes the components of the Gregorian Date (years, months, days, hours
    minutes, seconds) from a Julian Date. For reference, see Algorithm 22 in 
    Vallado (Fourth Edition) Section 3.6.6 pg 202
    
    Parameters
    ----------
    jd: double
        the Julian Date
    
    Returns
    -------
    (year, month, day, hours, minutes, seconds): tuple
        the year, month, day, hours, minutes, seconds of the date/time
    """
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