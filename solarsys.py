# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 11:27:13 2016

@author: Alex
"""

class Body:
    x = 0
    
class Moon(Body):
    
    # Orbit Params    
    sma_ER = 60.27
    sma_km = 384400.0
    ecc = 0.05490
    inc_deg = 5.145396
    raan_deg = None
    long_perihelion = None
    true_long = None
    period_yrs = 0.0748
    tropical_days = 27.321582
    orb_vel = 1.0232 #km/s
    
    # Body Params
    eq_radius_km = 1738.0
    flat = None
    mu = 4902.799
    mass_norm = 0.01230
    mass_kg = 7.3483E22
    rot_days = 27.32166
    eq_inc = 6.68
    j2 = 0.0002027
    j3 = None
    j4 = None
    density = 3.34
    
class Earth(Body):
    
    # Orbit Params
    sma_AU = 1.00000100178
    sma_km = 149598023.0
    ecc = 0.016708617
    inc = 0.0
    raan_deg = 0.0
    long_perihelion = 102.93734808
    true_long = 100.46644851
    period_yrs = 0.99997862
    
    # Body Params
    eq_radius_km = 6378.1363
    flag = 0.0033528131
    mu = 3.986004415E5
    mass_norm = 1.0
    mass_kg = 5.9742E24
    rot_days = 0.99726968
    eq_inc = 23.45
    j2 = 0.0010826269
    j3 = -0.0000025323
    j4 = -0.0000016204
    density = 5.515