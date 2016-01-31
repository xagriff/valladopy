# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 16:41:34 2016

@author: Alex
"""

import numpy as np

def find_c2c3(phi):
    if phi > 1E-6:
        sqrt_phi = np.sqrt(phi)
        sqrt_phi3 = np.sqrt(phi**3)
        c2 = (1 - np.cos(sqrt_phi))/phi
        c3 = (1 - np.sin(sqrt_phi))/sqrt_phi3
    else:
        if phi < -1E-6:
            sqrt_phi = np.sqrt(-phi)
            sqrt_phi3 = np.sqrt((-phi)**3)
            c2 = (1 - np.cosh(sqrt_phi))/phi
            c3 = (np.sinh(sqrt_phi) - sqrt_phi)/sqrt_phi3
        else:
            c2 = 0.5
            c3 = 1.0/6.0
    return c2, c3
     
def kep_eqtnE(M, e, tol=1E-8):
    if ((M > -np.pi) and (M < 0)) or M > np.pi:
        E = M - e
    else:
        E = M + e

    diff = np.inf
    while diff > tol:
        E_new = E + ((M - E + e*np.sin(E))/(1.0 - e*np.cos(E)))
        diff = np.abs(E_new - E)
        E = E_new
        
    return E
    


    