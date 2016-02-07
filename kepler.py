# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 16:41:34 2016

@author: Alex
"""

import numpy as np
from solarsys import *

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
        Mean Anomaly (radians)
    e: double
        Eccentricity
    tol: double, optional, default=1E-8
        Convergence tolerance used in Newton-Raphson method
        
    Returns
    -------
    E: double
        Eccentric Anomaly (radians)
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
    

def kep_eqtnP(del_t, p, mu=Earth.mu):
    """Parabolic solution to Kepler's Equation (Algorithm 3)
    
    A trigonometric approach to solving Kepler's Equation for 
    parabolic orbits. For reference, see Algorithm 3 in Vallado (Fourth 
    Edition), Section 2.2 (pg 69).
    
    Parameters
    ----------
    del_t: double
        Change in time
    p: double
        Semi-parameter
    mu: double, optional, default = Earth.mu
        Gravitational parameter; defaults to Earth
        
    Returns
    -------
    B: double
        Parabolic Anomaly (radians)
    """
    
    p3 = p**3
    n_p = 2.0*np.sqrt(mu/p3)
    s = 0.5*np.arctan(2.0/(3.0*n_p*del_t))
    w = np.arctan((np.tan(s))**(1/3.0))
    B = 2.0/np.tan(2.0*w)
    return B
    

def kep_eqtnH(M, e, tol=1E-8):
    """Hyperbolic solution to Kepler's Equation (Algorithm 4)
    
    A Newton-Raphson iterative approach to solving Kepler's Equation for 
    hyperbolic orbits. For reference, see Algorithm 4 in Vallado (Fourth 
    Edition), Section 2.2 (pg 71).
    
    Parameters
    ----------
    M: double
        Mean Anomaly (radians)
    e: double
        Eccentricity
    tol: double, optional, default=1E-8
        Convergence tolerance used in Newton-Raphson method
        
    Returns
    -------
    H: double
        Hyperbolic Anomaly (radians)
    """
    if e < 1.6:
        if (M > -np.pi and M < 0) or (M > np.pi):
            H = M-e
        else:
            H = M+e
    else:
        if e < 3.6 and (np.abs(M) > np.pi):
            H = M - np.sign(M)*e
        else:
            H = M/(e-1)
    
    diff = np.inf
    while diff > tol:
        H_new = H + (M - e*np.sinh(H) + H)/(e*np.cosh(H) -1)
        diff = np.abs(H_new - H)
        H = H_new
    return H
    

def true_to_anom(true_anom, e):
    """ Converts true anomaly to the proper orbit anomaly (Algorithm 5)
    
    Converts true anomaly to eccentric (E) anomaly for elliptical orbits, 
    parabolic (B) anomaly for parabolic orbits, or hyperbolic anomaly (H) for 
    hyperbolic orbits.  For reference, see Algorithm 5 in Vallado (Fourth 
    Edition), Section 2.2 (pg 77).
    
    Parameters
    ----------
    true_anom: double
        True anomaly (radians)
    e: double
        Eccentricity
        
    Returns
    -------
    E/B/H: double
        Eccentric, Parabolic, or Hyperbolic anomaly (radians)
    """
    if e < 1.0:
        num = np.sin(true_anom)*np.sqrt(1.0 - e**2)
        denom = 1 + e*np.cos(true_anom)
        E = np.arcsin(num/denom)
        return E
    elif e == 1.0:
        B = np.tan(0.5*true_anom)
        return B
    else:
        num = np.sin(true_anom)*np.sqrt(e**2 - 1)
        denom = 1 + e*np.cos(true_anom)
        H = np.arcsinh(num/denom)
        return H


def eccentric_to_true(E, e):
    """Converts Eccentric Anomaly (E) to true anomaly (Algorithm 6 (a))
    
    Converts eccentric anomaly (E) to true anomaly for elliptical orbits. For 
    reference, see Algorithm 6 in Vallado (Fourth Edition), Section 2.2 (pg 77)
    
    Parameters
    ----------
    E: double
        Eccentric anomaly (radians)
    e: double
        Eccentricity
    
    Returns
    -------
    true_anom: double
        True anomaly (radians)
    """
    num = np.cos(E) - e
    denom = 1 - e*np.cos(E)
    true_anom = np.arccos(num/denom)
    return true_anom
    

def parabolic_to_true(B, p, r):
    """Converts parabolic anomaly (B) to true anomaly (Algorithm 6 (b))
    
    Converts parabolic anomaly (B) to true anomaly for parabolic orbits. For
    reference, see Algorithm 6 in Vallado (Fourth Editon), Section 2.2 (pg 77)
    
    Parameters
    ----------
    B: double
        Parabolic anomaly (radians)
    p: double
        semi-parameter (km)
    r: double
        distance to the focus (km)
        
    Returns
    -------
    true_anom: double
        True anomaly (radians)
    """
    true_anom = np.arcsin(p*B/r)
    return true_anom
    

def hyperbolic_to_true(H, e):
    """Converts hyperbolic anomaly (H) to true anomaly (Algorithm 6 (c))
    
    Converts hyperbolic anomaly (H) to true anomaly for hyperbolic orbits. For
    reference, see Algorithm 6 in Vallado (Fourth Edition), Section 2.2 (pg 77)
    
    Parameters
    ----------
    H: double
        Hyperbolic anomaly (radians)
    e: double
        eccentricity
        
    Returns
    -------
    true_anom: double
        True anomaly (radians)
    """
    
    num = np.cosh(H) - e
    denom = 1 - e*np.cosh(H)
    true_anom = np.arccos(num/denom)
    return true_anom