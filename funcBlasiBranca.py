"""
Function based on Blasi and Branca 2001 kinetic reaction scheme for biomass
pyrolysis. Primary reactions evaluated at some temperature.

Reference:
Blasi and Branca 2001
"""

# Modules
# -----------------------------------------------------------------------------

import numpy as np

# Function
# -----------------------------------------------------------------------------

def blasibranca(T, pw, pg, pt, pc, dt, i):
    """
    INPUTS:
    T = temperature, K
    pw = wood concentration, kg/m^3
    pg = gas concentration, kg/m^3
    pt = tar concentration, kg/m^3
    pc = char concentration, kg/m^3
    dt = delta time, s
    i = index    
    OUTPUTS:
    pww = wood
    pgg = gas
    ptt = tar
    pcc = char
    """
    
    R = 0.008314 # universal gas constant, kJ/mol*K
    
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 4.38e9;    E1 = 152.7      # wood -> gas
    A2 = 1.08e10;   E2 = 148.0      # wood -> tar
    A3 = 3.27e6;    E3 = 111.7      # wood -> char
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char
    
    # reaction rate for each reaction, rho/s
    rww = -(K1+K2+K3) * pw[i-1]     # wood rate
    rwg = K1 * pw[i-1]              # wood -> gas rate
    rwt = K2 * pw[i-1]              # wood -> tar rate
    rwc = K3 * pw[i-1]              # wood -> char rate
    
    # wood, char, gas concentrations as a density, kg/m^3
    pww = pw[i-1] + rww*dt          # wood
    pgg = pg[i-1] + rwg*dt          # gas
    ptt = pt[i-1] + rwt*dt          # tar
    pcc = pc[i-1] + rwc*dt          # char
    
    # return the wood, char, gas, tar concentrations as a density, kg/m^3
    return pww, pgg, ptt, pcc
