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

# Function - primary kinetic reactions from Table 1
# -----------------------------------------------------------------------------

def blasi1(rhow, T, dt, nt):
    """
    rhow = wood density, kg/m^3
    T = temperature, K
    dt = time step, s
    nt = total number of time steps
    """
    
    # vector for initial wood concentration, kg/m^3
    pw = np.ones(nt)*rhow
    
    # vectors to store product concentrations, kg/m^3
    pg = np.zeros(nt)    # gas
    pt = np.zeros(nt)    # tar
    pc = np.zeros(nt)    # char
    
    R = 0.008314 # universal gas constant, kJ/mol*K
    
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 5.16e6;     E1 = 88.6     # wood -> gas
    A2 = 1.48e10;    E2 = 112.7    # wood -> tar
    A3 = 2.66e10;    E3 = 106.5    # wood -> char
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char
    
    # concentrations at each time step for each product, kg/m^3
    # reaction rate as r, rho/s
    # concentration as density p, kg/m^3
    for i in range(1, nt):
        rww = -(K1+K2+K3) * pw[i-1]     # wood rate
        rwg = K1 * pw[i-1]              # wood -> gas rate
        rwt = K2 * pw[i-1]              # wood -> tar rate
        rwc = K3 * pw[i-1]              # wood -> char rate
        pw[i] = pw[i-1] + rww*dt        # wood
        pg[i] = pg[i-1] + rwg*dt        # gas
        pt[i] = pt[i-1] + rwt*dt        # tar
        pc[i] = pc[i-1] + rwc*dt        # char
    
    # return the wood, gas, tar, char concentrations as a density, kg/m^3
    return pw, pg, pt, pc


# Function - primary and secondary reactions from Table 1
# -----------------------------------------------------------------------------

def blasi2(rhow, T, dt, nt):
    """
    rhow = wood density, kg/m^3
    T = temperature, K
    dt = time step, s
    nt = total number of time steps
    """
    
    # vector for initial wood concentration, kg/m^3
    pw = np.ones(nt)*rhow
    
    # vectors to store product concentrations, kg/m^3
    pg = np.zeros(nt)    # gas
    pt = np.zeros(nt)    # tar
    pc = np.zeros(nt)    # char
    
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
    
    # concentrations at each time step for each product, kg/m^3
    # reaction rate as r, rho/s
    # concentration as density p, kg/m^3
    for i in range(1, nt):
        rww = -(K1+K2+K3) * pw[i-1]     # wood rate
        rwg = K1 * pw[i-1]              # wood -> gas rate
        rwt = K2 * pw[i-1]              # wood -> tar rate
        rwc = K3 * pw[i-1]              # wood -> char rate
        rtg = K4 * pt[i-1]              # tar -> gas
        rtc = K5 * pt[i-1]              # tar -> char
        pw[i] = pw[i-1] + rww*dt              # wood
        pg[i] = pg[i-1] + (rwg+rtg)*dt        # gas
        pt[i] = pt[i-1] + (rwt-rtg-rtc)*dt    # tar
        pc[i] = pc[i-1] + (rwc+rtc)*dt        # char
    
    # return the wood, gas, tar, char concentrations as a density, kg/m^3
    return pw, pg, pt, pc
    