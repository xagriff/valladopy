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
        
        
        
def test_main():
        support.run_unittest(TimeExamplesFromBookTestCase)
        
if __name__ == '__main__':
        test_main()