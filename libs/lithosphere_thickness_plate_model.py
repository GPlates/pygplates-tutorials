import numpy as np
import math

def plate_isotherm_depth(age, temp, *vartuple) :
    "Computes the depth to the temp - isotherm in a cooling plate mode.\
    Solution by iteration. By default the plate thickness is 125 km as\
    in Parsons/Sclater.  Change given a 3rd parameter."

    if len(vartuple) != 0 :
        PLATE_THICKNESS_KM = vartuple[0]
    else :
        PLATE_THICKNESS_KM = 125

    PLATE_THICKNESS = PLATE_THICKNESS_KM * 1000
    
    
    # default depth is 0
    z = 0

    if age <= 0.0 :
        z_try = 0
        done = 1
    else :
        z_too_small = 0.0
        z_too_big = PLATE_THICKNESS
        done = 0
        n_try = 0
            
    while done != 1 and n_try < 20 :
        n_try += 1
        z_try = 0.5 * (z_too_small + z_too_big)
        t_try = plate_temp (age, z_try, PLATE_THICKNESS)
        t_wrong = temp - t_try

        if t_wrong < -0.001 :
            z_too_big = z_try
        elif t_wrong > 0.001 :
            z_too_small = z_try
        else :
            done = 1

        z = z_try
    return z

def plate_temp(age, z, PLATE_THICKNESS) :
    "Computes the temperature in a cooling plate for age = t\
    and at a depth = z."

    KAPPA = 0.804E-6
    T_MANTLE = 1350.0
    T_SURFACE = 0.0
    SEC_PR_MY = 3.15576e13

    t = T_SURFACE

    sum = 0
    sine_arg = math.pi * z / PLATE_THICKNESS
    exp_arg = -KAPPA * math.pi * math.pi * age * SEC_PR_MY / (PLATE_THICKNESS * PLATE_THICKNESS)
    for k in range(1, 20) :
        sum = sum + np.sin(k * sine_arg) * np.exp(k*k*exp_arg)/k

    if age <= 0.0 :
        t = T_MANTLE * np.ones(z.shape)
    else :
        t = t + 2.0 * sum * (T_MANTLE - T_SURFACE)/math.pi + (T_MANTLE - T_SURFACE) * z/PLATE_THICKNESS
    
    return t