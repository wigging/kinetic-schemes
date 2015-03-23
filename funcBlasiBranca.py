"""
Function based on Blasi 2001 kinetic reaction scheme for biomass pyrolysis. 
Reactions are evaluated at some temperature.

Reference:
Blasi, Branca, 2001. Ind. Eng. Chem. Res., 40, pp.5547-5556.
"""

# Modules
# -----------------------------------------------------------------------------

import numpy as np

# Function - primary kinetic reactions from Table 1
# -----------------------------------------------------------------------------

def blasibranca(rhow, T, dt, nt):
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
    A1 = 4.38e9;    E1 = 152.7    # wood -> gas
    A2 = 1.08e10;   E2 = 148      # wood -> tar
    A3 = 3.27e6;    E3 = 111.7    # wood -> char
    
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
        pw[i] = pw[i-1] + rww*dt          # wood
        pg[i] = pg[i-1] + rwg*dt          # gas
        pt[i] = pt[i-1] + rwt*dt          # tar
        pc[i] = pc[i-1] + rwc*dt          # char 
    
    # return the wood, gas, tar, char concentrations as a density, kg/m^3
    return pw, pg, pt, pc
    
    