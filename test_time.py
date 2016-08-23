# -*- coding: utf-8 -*-
"""
Created on Sat Aug 20 17:22:44 2016

@author: Alex
"""

import unittest
from test import support
import astro_time
import math
import numpy as np

class TimeExamplesFromBookTestCase(unittest.TestCase):
    
    def test_example_3_4_julian_date(self):
        month = 10
        day = 26
        yr = 1996
        hr = 14
        minute = 20
        sec = 0.0
        act_jd = astro_time.julian_date(yr, month, day, hr, minute, sec)
        exp_jd = 2450383.09722222
        self.assertAlmostEqual(act_jd, exp_jd, 8)
        
    def test_example_3_5_gmst_lst(self):
        month = 8
        day = 20
        yr = 1992
        hr = 12
        minute = 14
        sec = 0.0
        lon = -104.0
        
        jd_ut1 = astro_time.julian_date(yr, month, day, hr, minute, sec)
        act_theta_gmst = astro_time.find_gmst(jd_ut1)
        act_theta_lst = astro_time.find_lst(act_theta_gmst, lon)
        
        exp_theta_gmst = 152.578787810
        exp_theta_lst = 48.578787810
        
        self.assertAlmostEqual(act_theta_gmst, exp_theta_gmst)
        self.assertAlmostEqual(act_theta_lst, exp_theta_lst)
        
    def test_example_3_8_dms2rad(self):
        degrees = -35.0
        minutes = -15.0
        seconds = -53.63
        
        act_rad = astro_time.dms2rad(degrees, minutes, seconds)
        exp_rad = -0.6154886
        
        self.assertAlmostEqual(act_rad, exp_rad)
        
    def test_example_3_8_rad2dms(self):
        rad = -0.6154886
        
        (act_deg, act_min, act_sec) = astro_time.rad2dms(rad)
        exp_deg = -35.0
        exp_min = -15.0
        exp_sec = -53.6368
        
        self.assertEqual(act_deg, exp_deg)
        self.assertEqual(act_min, exp_min)
        self.assertAlmostEqual(act_sec, exp_sec, 4)
        
    def test_example_3_9_hms2rad(self):
        hours = 15
        minutes = 15
        seconds = 53.63
        
        act_rad = astro_time.hms2rad(hours, minutes, seconds)
        exp_rad = 3.996341
        
        self.assertAlmostEqual(act_rad, exp_rad, 6)
        
    def test_example_3_9_rad2hms(self):
        rad = 3.996341
        (act_hrs, act_min, act_sec) = astro_time.rad2hms(rad)
        exp_hrs = 15
        exp_min = 15
        exp_sec = 53.63
        
        self.assertEqual(act_hrs, exp_hrs)
        self.assertEqual(act_min, exp_min)
        self.assertAlmostEqual(act_sec, exp_sec, 2)
        
    def test_example_3_10_time2hms(self):
        tau = 48165.98
        (act_hrs, act_min, act_sec) = astro_time.time2hms(tau)
        exp_hrs = 13
        exp_min = 22
        exp_sec = 45.98
        
        self.assertEqual(act_hrs, exp_hrs)
        self.assertEqual(act_min, exp_min)
        self.assertAlmostEqual(act_sec, exp_sec, 2)
        
    def test_example_3_10_hms2time(self):
        hours = 13
        minutes = 22
        seconds = 45.98
        act_time = astro_time.hms2time(hours, minutes, seconds)
        exp_time = 48165.98
        
        self.assertAlmostEqual(act_time, exp_time, 2)
        
class OtherTimeFunctionsTestCase(unittest.TestCase):
    
    def test_is_leap_year(self):
        
        y1 = 1996
        y2 = 1900
        y3 = 2000
        y4 = 2005
        y5 = 2100
        
        b1 = astro_time.is_leap_year(y1)
        b2 = astro_time.is_leap_year(y2)
        b3 = astro_time.is_leap_year(y3)
        b4 = astro_time.is_leap_year(y4)
        b5 = astro_time.is_leap_year(y5)
        
        self.assertEqual(b1, True)
        self.assertEqual(b2, False)
        self.assertEqual(b3, True)
        self.assertEqual(b4, False)
        self.assertEqual(b5, False)
        
        
def test_main():
        support.run_unittest(TimeExamplesFromBookTestCase)
        support.run_unittest(OtherTimeFunctionsTestCase)
        
if __name__ == '__main__':
        test_main()