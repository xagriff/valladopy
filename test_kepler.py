# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 22:50:18 2016

@author: Alex
"""

import unittest
from test import support
import kepler
import math
import numpy as np

class ExamplesFromBookTestCase(unittest.TestCase):
    
#    def create_assertEqual_msg(self, actual, expected):
#        message = ('The expected value was {expected_str}, '
#            'but the function returned {actual_str}').format(
#            expected_str=expected, actual_str=actual)
#        return message
    
    def test_example_2_1_Keplers_Equation(self):
        M = 235.4
        e = 0.4
        M_rad = math.radians(M)
        actual = math.degrees(kepler.kep_eqtnE(M_rad, e))
        expected = 220.512074767522
        self.assertAlmostEqual(actual, expected, places=8)
        
    def test_example_2_2_Keplers_Equation_Parabolic(self):
        del_t = 53.7874
        p = 25512.0
        
        del_t_sec = del_t*60.0
        actual = kepler.kep_eqtnP(del_t_sec, p)
        expected = 0.817751
        self.assertAlmostEqual(actual, expected, places=6)
        
    def test_example_2_3_Keplers_Equation_Hyperbolic(self):
        M = 235.4
        e = 2.4
        M_rad = math.radians(M)
        actual = kepler.kep_eqtnH(M_rad, e)
        expected = 1.601376144
        self.assertAlmostEqual(actual, expected, places=8)
        
    def test_example_2_4_Keplers_Problem(self):
        r_ijk = np.matrix((1131.340, -2282.343, 6672.423)).T
        v_ijk = np.matrix((-5.64305, 4.30333, 2.42879)).T
        del_t = 40.0
        
        del_t_sec = del_t*60.0
        
        (actual_pos, actual_vel) = kepler.kepler(r_ijk, v_ijk, del_t_sec)
        expected_pos = np.matrix((-4219.7527, 4363.0292, -3958.7666)).T
        expected_vel = np.matrix((3.689866, -1.916735, -6.112511)).T
        self.assertAlmostEqual(actual_pos.item(0), expected_pos.item(0), places=3)
        self.assertAlmostEqual(actual_pos.item(1), expected_pos.item(1), places=3)
        self.assertAlmostEqual(actual_pos.item(2), expected_pos.item(2), places=3)
        self.assertAlmostEqual(actual_vel.item(0), expected_vel.item(0), places=5)
        self.assertAlmostEqual(actual_vel.item(1), expected_vel.item(1), places=5)
        self.assertAlmostEqual(actual_vel.item(2), expected_vel.item(2), places=5)
                               
    def test_example_2_5_find_orbital_elements(self):
        r_ijk = np.matrix((6524.834, 6862.875, 6448.296)).T
        v_ijk = np.matrix((4.901327, 5.533756, -1.976341)).T
        (pe, ae, ee, ie, raane, aope, tanome) = (11067.790, 36127.343, 0.832853, math.radians(87.870), math.radians(227.898), math.radians(53.38), math.radians(92.335))
        (pa, aa, ea, ia, raana, aopa, tanoma) = kepler.rv2coe(r_ijk, v_ijk)
        self.assertAlmostEqual(pa, pe, places=1)
        self.assertAlmostEqual(aa, ae, places=1)
        self.assertAlmostEqual(ea, ee, places=6)
        self.assertAlmostEqual(ia, ie, places=3)
        self.assertAlmostEqual(raana, raane,places=3)
        self.assertAlmostEqual(aopa, aope, places=3)
        self.assertAlmostEqual(tanoma, tanome, places=3)
         
def test_main():
        support.run_unittest(ExamplesFromBookTestCase)
        
if __name__ == '__main__':
        test_main()