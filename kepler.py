# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 16:41:34 2016

@author: Alex
"""

import numpy as np

def find_c2c3(psi):
    """c(:math:`\psi`) functions for the universal formulation (Algorithm 1)
    
    A trigonometric implementation of the :math:`c(\psi)` functions needed in 
    the universal formulation of Kepler's Equation. For reference, see 
    Algorithm 1 in Vallado (Fourth Edition), Section 2.2 (pg 63).
    
    Parameters
    ----------
    psi: double
        :math:`\psi = \chi^2/a` where :math:`\chi` is universal 
        variable and a is the semi-major axis
        
    Returns
    -------
    c2: double
        c2 coefficient in universal formulation of Kepler's Equation
    c3: double, coefficient
        c3 coefficient in universal formulation of Kepler's Equation
    """

    if psi > 1E-6:
        sqrt_psi = np.sqrt(psi)
        sqrt_psi3 = np.sqrt(psi**3)
        c2 = (1 - np.cos(sqrt_psi))/psi
        c3 = (1 - np.sin(sqrt_psi))/sqrt_psi3
    else:
        if psi < -1E-6:
            sqrt_psi = np.sqrt(-psi)
            sqrt_psi3 = np.sqrt((-psi)**3)
            c2 = (1 - np.cosh(sqrt_psi))/psi
            c3 = (np.sinh(sqrt_psi) - sqrt_psi)/sqrt_psi3
        else:
            c2 = 0.5
            c3 = 1.0/6.0
    return c2, c3
     
def kep_eqtnE(M, e, tol=1E-8):
    """Elliptical solution to Kepler's Equation (Algorithm 2)
    
    A Newton-Raphson iterative approach to solving Kepler's Equation for 
    elliptical orbits. For reference, see Algorithm 2 in Vallado (Fourth 
    Edition), Section 2.2 (pg 65).
    
    Parameters
    ----------
    M: double
        Mean Anomaly
    e: double
        Eccentricity
    tol: double, optional, default=1E-8
        Convergence tolerance used in Newton-Raphson method
        
    Returns
    -------
    E: double
        Eccentric Anomaly
    """
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
    


    