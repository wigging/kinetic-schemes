"""
Function based on Janse 2000 kinetic reaction scheme for biomass pyrolysis. 
Primary and secondary reactions evaluated at some temperature.

Functions:
janse1 - primary reactions only
janse2 - primary and secondary reactions

Reference:
Janse, Westerhout, Prins, 2000. Chem. Eng. Process., 39, 239-252.
"""

# Modules
# -----------------------------------------------------------------------------

import numpy as np

# Function - primary reactions from Wagenaar 1993 (1-3), see Table 1
# -----------------------------------------------------------------------------

def janse1(T, pw, pg, pt, pc, dt, i):
    """
    Primary kinetic reactions from Wagenaar 1993 (1-3), see Table 1.
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
    A1 = 1.11e11;   E1 = 177    # wood -> gas
    A2 = 9.28e9;    E2 = 149    # wood -> tar
    A3 = 3.05e7;    E3 = 125    # wood -> char
    
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
    
    # return the wood, gas, tar, char concentrations as a density, kg/m^3
    return pww, pgg, ptt, pcc
    
    
# Function - primary and secondary reactions from Liden 1988 (4-5), see Table 1
# -----------------------------------------------------------------------------

def janse2(T, pw, pg, pt, pc, dt, i):
    """
    Primary and secondary reactions from Liden 1998 (4-5), see Table 1.
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
    A1 = 1.11e11;   E1 = 177     # wood -> gas
    A2 = 9.28e9;    E2 = 149     # wood -> tar
    A3 = 3.05e7;    E3 = 125     # wood -> char
    A4 = 8.6e4;     E4 = 87.8    # tar -> gas
    A5 = 7.7e4;     E5 = 87.8    # tar -> char

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
    rtg = K4 * pt[i-1]              # tar -> gas rate
    rtc = K5 * pt[i-1]              # tar -> char rate
    
    # wood, char, gas concentrations as a density, kg/m^3
    pww = pw[i-1] + rww*dt                  # wood
    pgg = pg[i-1] + (rwg + rtg)*dt          # gas
    ptt = pt[i-1] + (rwt - rtg - rtc)*dt    # tar
    pcc = pc[i-1] + (rwc + rtc)*dt          # char
    
    # return the wood, gas, tar, char concentrations as a density, kg/m^3
    return pww, pgg, ptt, pcc
    
    