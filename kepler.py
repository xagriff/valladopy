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
        psi3 = psi*psi*psi
        sqrt_psi3 = np.sqrt(psi3)
        c2 = (1.0 - np.cos(sqrt_psi))/psi
        c3 = (sqrt_psi - np.sin(sqrt_psi))/sqrt_psi3
    else:
        if psi < -1E-6:
            sqrt_psi = np.sqrt(-psi)
            sqrt_psi3 = np.sqrt((-psi)*(-psi)*(-psi))
            c2 = (1.0 - np.cosh(sqrt_psi))/psi
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
    

def parabolic_to_true(B):#, p, r):
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
    #true_anom = np.arcsin(p*B/r)
    true_anom = np.arctan(B)
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
    
def rv2coe(r, v, mu=Earth.mu):
    """Converts position/velocity to Keplerian orbital elements (Algorithm 9)
    
    Converts position and velocity vectors in the IJK frame to Keplerian 
    orbital elements.  For reference, see Algorithm 9 in Vallado (Fourth 
    Edition), Section 2.5 (pg 113)
    
    Parameters
    ----------
    r: numpy.matrix (3x1)
        Position vector (km)
    v: numpy.matrix (3x1)
        Velocity vector (km/s)
    mu: double, optional, default=3.986004415E5 (Earth.mu in solarsys.py)
        Gravitational parameter (km^3/s^2)
        
    Returns
    -------
    p: double
        Semi-parameter (km)
    a: double
        Semi-major axis (km)
    ecc: double
        Eccentricity
    inc: double
        Inclination (radians)
    raan: double
        Right ascension of the ascending node (radians)
    aop: double
        Argument of perigee (radians)
    t_anom: double
        True anomaly (radians)
    special: string
        String indicating special case (circular, equatorial, etc.). If not
        special, returned as None
        
    Note
    ----
    This algorithm handles special cases (circular, equatorial, etc.) by 
    setting raan, aop, and anom as would be used by coe2rv (Algorithm 10).
    """
    r_mag = np.linalg.norm(r)
    r_inv = 1.0/r_mag # Store for efficiency
    v_mag = np.linalg.norm(v)
    h = np.matrix(np.cross(r, v, axis=0)) # Column vectors in, column vec out
    h_mag = np.linalg.norm(h)
    k_hat = np.matrix([[0.0],[0.0],[1.0]])
    n = np.matrix(np.cross(k_hat, h, axis=0))
    n_mag = np.linalg.norm(n)
    e_vec = ((v_mag*v_mag - mu*r_inv)*r - ((r.T*v)[0,0]*v))/mu
    ecc = np.linalg.norm(e_vec)
    e_inv = 1.0/ecc # Store for efficiency
    xi = 0.5*v_mag*v_mag - mu*r_inv
    if ecc != 1.0:
        a = -0.5*mu/xi
        p = a*(1.0 - ecc*ecc)
    else:
        a = np.inf
        p = h_mag*h_mag/mu
    inc = np.arccos(h[2,0]/h_mag)
    
    if n_mag == 0.0: # Equatorial
        if ecc == 0.0: # Circular
            lambda_true = np.arccos(r[0,0]*r_inv)
            if r[1] < 0:
                lambda_true = 2.0*np.pi - lambda_true
            raan = 0.0
            aop = 0.0
            t_anom = lambda_true
            return p, a, ecc, inc, raan, aop, t_anom
        else:
            e_inv = 1.0/ecc
            omega_true = np.arccos(e_vec[0,0]*e_inv)
            if e_vec[1] < 0:
                omega_true = 2.0*np.pi - omega_true
            raan = 0.0
            aop = omega_true
            t_anom = np.arccos(np.dot(e_vec, r)*e_inv*r_inv)
            return p, a, ecc, inc, raan, aop, t_anom
    elif ecc == 0.0: # Circular
        n_inv = 1.0/n_mag
        raan = np.arccos(n[0,0]*n_inv)
        aop = 0.0
        u = np.arccos(np.dot(n, r)*n_inv*r_inv)
        if r[2] < 0:
            u = 2.0*np.pi - u
        t_anom = u
        return p, a, ecc, inc, raan, aop, t_anom
    else:
        n_inv = 1.0/n_mag
        e_inv = 1.0/ecc
        raan = np.arccos(n[0,0]*n_inv)
        if n[1] < 0:
            raan = 2.0*np.pi - raan
        aop = np.arccos((n.T*e_vec)[0,0]*n_inv*e_inv)
        if e_vec[2] < 0:
            aop = 2.0*np.pi - aop
        t_anom = np.arccos((e_vec.T*r)[0,0]*e_inv*r_inv)
        if (r.T*v)[0,0] < 0:
            t_anom = 2.0*np.pi - t_anom
        return p, a, ecc, inc, raan, aop, t_anom
    
def coe2rv(p, ecc, inc, raan, aop, anom, mu=Earth.mu):
    """Converts Keplerian orbital elements to pos/vel vectors (Algorithm 10)
    
    Converts Keplerian orbital elements to position/velocity vectors (km, km/s)
    in the IJK frame.  For reference, see Algorithm 10 in Vallado (Fourth 
    Edition), Section 2.6 (pg 118).
    
    Parameters
    ----------
    p: double
        Semi-parameter (km)
    ecc: double
        Eccentricity
    inc: double
        Inclination (radians)
    raan: double
        Right Ascension of the Ascending Node (radians)
    aop: double
        Argument of perigee (radians)
    anom: double
        True anomaly (radians)
    mu: double, optional, default=3.986004415E5 (Earth.mu in solarsys.py)
        Gravitational parameter (km^3/s^2)
    
    Returns
    -------
    r_ijk: numpy.matrix (3x1)
        Position vector in the IJK frame (km)
    v_ijk: numpy.matrix (3x1)
        Velocity vector in the IJK frame (km/s)
    
    Note
    ----
    Algorithm assumes that raan, aop, and anom have been set to account for 
    special cases (circular, equatorial, etc.) as in rv2coe (Algorithm 9)
    """
    # Stored trig comps
    cosv = np.cos(anom)
    sinv = np.sin(anom)
    cosi = np.cos(inc)
    sini = np.sin(inc)
    cosw = np.cos(aop)
    sinw = np.sin(aop)
    coso = np.cos(raan)
    sino = np.sin(raan)
    
    r_pqw = np.matrix([p*cosv/(1.0+ecc*cosv), p*sinv/(1.0+ecc*cosv), 0.0])
    r_pqw = r_pqw.T # Make column vector
    v_pqw = np.matrix([-np.sqrt(mu/p)*sinv, np.sqrt(mu/p)*(ecc+cosv), 0.0])
    v_pqw = v_pqw.T # Make column vector
    
    m_pqw2ijk = [[coso*cosw-sino*sinw*cosi, 
                             -coso*sinw-sino*cosw*cosi, sino*sini],
                        [sino*cosw+coso*sinw*cosi, -sino*sinw+coso*cosw*cosi, 
                         -coso*sini], [sinw*sini, cosw*sini, cosi]]
    m_pqw2ijk = np.matrix(m_pqw2ijk)
#    m_pqw2ijk = np.matrix([[row1], [row2], [row3]])
    # Convert to IJK frame
    r_ijk = m_pqw2ijk*r_pqw
    v_ijk = m_pqw2ijk*v_pqw
    
    return r_ijk, v_ijk
    
def findTOF(r0, r, p, mu=Earth.mu):
    """Finds the Time of Flight between two position vectors (Algorithm 11)
    
    Finds the time of flight bewteen two position vectors in the IJK frame.  
    For reference, see Algorithm 11 in Vallado (Fourth Edition), Section 2.8
    (pg 126)
    
    Parameters
    ----------
    r0: numpy.matrix (3x1)
        Initial position vector (km)
    r: numpy.matrix (3x1)
        Second position vector (km)
    p: double
        Semi-parameter (km)
    mu: double, optional, default=3.986004415E5 (Earth.mu in solarsys.py)
        Gravitational parameter (km^3/s^2)
        
    Returns
    -------
    TOF: double
        Time of flight (seconds)
    """

    r0_mag = np.linalg.norm(r0)
    r_mag = np.linalg.norm(r)
    
    cosdv = np.dot(r0.T, r)/(r0_mag*r_mag)
    del_anom = np.arccos(cosdv)
    k = r0_mag*r_mag*(1.0 - cosdv)
    l = r0_mag + r_mag
    m = r0_mag*r_mag*(1.0 + cosdv)
    
    a = m*k*p/((2.0*m - l*l)*p*p + 2.0*k*l*p - k*k)
    f = 1.0 - (r_mag/p)*(1.0 - cosdv)
    g = r0_mag*r_mag*np.sin(del_anom)/(np.sqrt(mu*p))
    
    if a > 0.0:
        if a==np.inf:
            c = np.sqrt(r0_mag**2 + r_mag**2 - 2.0*r0_mag*r_mag*cosdv)
            s = (r0_mag + r_mag + c)*0.5
            TOF = (2.0/3.0)*np.sqrt(0.5*s**3/mu)*(1.0 - ((s - c)/s)**1.5)
            return TOF
        else:
            f_dot = np.sqrt(mu/p)*np.tan(0.5*del_anom)*((1.0 - cosdv)/p - 
                1.0/r0_mag - 1/r_mag)
            cosde = 1.0 - (r0_mag/a)*(1.0 - f)
            sinde = (-r0_mag*r_mag*f_dot)/(np.sqrt(mu*a))
            del_E = np.arccos(cosde)
            TOF = g + np.sqrt(a**3/mu)*(del_E - sinde)
            return TOF
    elif a < 0.0:
        coshdh = 1.0 + (f - 1.0)*(r0_mag/a)
        del_H = np.arccosh(coshdh)
        TOF = g + np.sqrt((-a)**3/mu)*(np.sinh(del_H) - del_H)
        return TOF
    else:
        TOF = None
        return TOF
        
def keplerCOE(r0, v0, dt, mu=Earth.mu):
    """Two body orbit propagation using classical orbital elements (Algorithm 7)
    
    Two body orbit propagation that uses a change to classical orbital elements
    and an update to true anomaly to find the new position and velocity 
    vectors. For reference, see Algorithm 7 in Vallado Section 2.3.1 (pg 81).
    
    Parameters
    ----------
    r0: list, numpy.array, or numpyt.matrix (length 3)
        Initial position vector (km)
    v0: list, numpy.array, or numpyt.matrix (length 3)
        Initial velocity vector (km/s)
    dt: double
        Time to propagate (seconds)
    mu: double, optional, default=3.986004415E5 (Earth.mu in solarsys.py)
        Gravitational parameter (km^3/s^2)
    
    Returns
    -------
    r: numpy.matrix (3x1)
        New position vector (km)
    v: numpy.matrix (3x1)
        New velocity vector (km/s)    
    """
    [p, a, ecc, inc, raan, aop, t_anom] = rv2coe(r0, v0, mu)
    n = np.sqrt(mu/a**3)
    
    if ecc != 0.0:
        anom0 = true_to_anom(t_anom, ecc)
    else:
        anom0 = t_anom
    
    if ecc < 1.0: # ecc < 1 -> elliptical or circular
        M0 = anom0 - ecc*np.sin(anom0)
        M = M0 + n*dt
        E = kep_eqtnE(M, ecc)
        if ecc != 0.0:
            t_anom = eccentric_to_true(E, ecc)
        else: # circular
            t_anom = E
    elif e == 1.0: # ecc = 1.0 -> Parabolic
        M0 = anom0 + anom0**3/3
        B = kep_eqtnP(dt, p)
        t_anom = parabolic_to_true(B)
    else: # ecc > 1.0 -> hyperbolic
        M0 = ecc*np.sinh(anom0) - anom0
        M = M0 + n*dt
        H = kep_eqtnH(M, ecc)
        t_anom = hyperbolic_to_true(H, ecc)
        
    r, v = coe2rv(p, ecc, inc, raan, aop, t_anom)
    return r, v
        
def kepler(r0, v0, dt, mu=Earth.mu, tol=1E-6):
    """Two body orbit propagation using universal variables (Algorithm 8)
    
    Two body orbit propagation that uses the universal formulation to find the 
    new position and velocity vectors. For reference, see Algorithm 8 in 
    Vallado Section 2.3.1 (pg 93).
    
    Parameters
    ----------
    r0: numpy.matrix (3x1)
        Initial position vector (km)
    v0: numpy.matrix (3x1)
        Initial velocity vector (km/s)
    dt: double
        Time to propagate (seconds)
    mu: double, optional, default=3.986004415E5 (Earth.mu in solarsys.py)
        Gravitational parameter (km^3/s^2)
    tol: double, optional, default=1E-6
        Convergence criterion (iteration difference)
    
    Returns
    -------
    r: numpy.matrix (3x1)
        New position vector (km)
    v: numpy.matrix (3x1)
        New velocity vector (km/s)    
    """
    
    r_mag = np.linalg.norm(r0)
    v_mag = np.linalg.norm(v0)
    alpha = -v_mag**2/mu + 2.0/r_mag
    
    if alpha > 0.000001: # ELliptical or circular
        chi0 = np.sqrt(mu)*dt*alpha
    elif alpha < -0.000001: # Hyperbolic
        a = 1/alpha
        chi0 = np.sign(dt)*np.sqrt(-a)*np.log(-2.0*mu*alpha*dt/(np.dot(r0.T,v0)
            + np.sign(dt)*np.sqrt(-mu*a)*(a - r_mag*alpha)))
    else: # Parabolic
        h = np.cross(r0, v0, axis=0)
        h_mag = np.linalg.norm(h)
        p = h_mag*h_mag/mu
        s = 0.5*np.arctan(1.0/(3.0*np.sqrt(mu/p)*dt))
        w = np.arctan((np.tan(s))**(1/3.0))
        chi0 = np.sqrt(p)*2.0/np.tan(2.0*w)
        
    diff = np.inf
    chi = chi0
    
    while np.abs(diff) > tol:
        psi = chi*chi*alpha
        c2, c3 = find_c2c3(psi)
        r = (chi**2*c2 + (np.dot(r0.T,v0)/np.sqrt(mu))*chi*(1.0 - psi*c3) + 
            r_mag*(1.0 - psi*c2))
        chi_new = chi + (np.sqrt(mu)*dt - chi**3*c3 - 
            (np.dot(r0.T, v0)/np.sqrt(mu))*chi**2*c2 - 
            r_mag*chi*(1.0 - psi*c3))/r
        diff = np.abs(chi_new - chi)
        chi = chi_new
    
    f = 1.0 - (chi**2/r_mag)*c2
    g = dt - (chi**3/np.sqrt(mu))*c3
    g_dot = 1.0 - (chi**2/r)*c2
    f_dot = (np.sqrt(mu)/(r*r_mag))*chi*(psi*c3 - 1.0)
    r_vec = f*r0 + g*v0
    v_vec = f_dot*r0 + g_dot*v0
    
    return r_vec, v_vec
   
