"""
Functions based on Blasi 1993 kinetic reaction scheme for biomass pyrolysis. 
Reactions are evaluated at some temperature.

Functions:
blasi1 - primary reactions only
blasi2 - primary and secondary reactions

Reference:
Blasi, 1993. Combustion Science and Technology, 90, pp.315–340.
"""

# Modules
# -----------------------------------------------------------------------------

import numpy as np

# Function
# -----------------------------------------------------------------------------

def blasi1(T, pw, pg, pt, pc, dt, i):
    """
    Primary kinetic reactions from Table 1.
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
    A1 = 5.16e6;     E1 = 88.6     # wood -> gas
    A2 = 1.48e10;    E2 = 112.7    # wood -> tar
    A3 = 2.66e10;    E3 = 106.5    # wood -> char
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char
    
    # reaction rate for each reaction, rho/s
    rww = -(K1+K2+K3) * pw[i-1]     # wood rate
    rwg = K1 * pw[i-1]              # wood -> gas rate
    rwt = K2 * pw[i-1]              # wood -> tar rate
    rwc = K3 * pw[i-1]              # wood -> char rate
    
    # wood, gas, tar, char concentrations as a density, kg/m^3
    pww = pw[i-1] + rww*dt          # wood
    pgg = pg[i-1] + rwg*dt          # gas
    ptt = pt[i-1] + rwt*dt          # tar
    pcc = pc[i-1] + rwc*dt          # char
    
    # return the wood, gas, tar, char concentrations as a density, kg/m^3
    return pww, pgg, ptt, pcc


# Function
# -----------------------------------------------------------------------------

def blasi2(T, pw, pg, pt, pc, dt, i):
    """
    Primary and secondary reations from Table 1.
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
    A1 = 5.16e6;     E1 = 88.6     # wood -> gas
    A2 = 1.48e10;    E2 = 112.7    # wood -> tar
    A3 = 2.66e10;    E3 = 106.5    # wood -> char
    A4 = 4.28e6;     E4 = 108      # tar -> gas
    A5 = 1.0e6;      E5 = 108      # tar -> char

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas
    K5 = A5 * np.exp(-E5 / (R * T))  # tar -> char

    # reaction rate for each reaction, rho/s
    rww = -(K1+K2+K3) * pw[i-1]     # wood rate
    rwg = K1 * pw[i-1]              # wood -> gas rate
    rwt = K2 * pw[i-1]              # wood -> tar rate
    rwc = K3 * pw[i-1]              # wood -> char rate
    rtg = K4 * pt[i-1]              # tar -> gas
    rtc = K5 * pt[i-1]              # tar -> char
    
    # wood, gas, tar, char concentrations as a density, kg/m^3
    pww = pw[i-1] + rww*dt              # wood
    pgg = pg[i-1] + (rwg+rtg)*dt        # gas
    ptt = pt[i-1] + (rwt-rtg-rtc)*dt    # tar
    pcc = pc[i-1] + (rwc+rtc)*dt        # char
    
    # return the wood, gas, tar, char concentrations as a density, kg/m^3
    return pww, pgg, ptt, pcc
    