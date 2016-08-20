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
        print(act_jd)
        exp_jd = 2450383.09722222
        self.assertAlmostEqual(act_jd, exp_jd, 8)
        
        
def test_main():
        support.run_unittest(TimeExamplesFromBookTestCase)
        
if __name__ == '__main__':
        test_main()