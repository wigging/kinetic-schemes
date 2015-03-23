"""
Function based on Font 1990 kinetic reaction scheme for biomass pyrolysis. 
Reactions evaluated at some temperature.

Functions:
font1 - fluidized bed kinetics
font2 - pyroprobe kinetics

Reference:
Font, Marcilla, Verdu, Devesa, 1990. Ind. Eng. Chem. Res., 29, pp.1846-1855.
"""

# Modules
# -----------------------------------------------------------------------------

import numpy as np

# Function - primary reactions from fluidized bed
# -----------------------------------------------------------------------------

def font1(rhow, T, dt, nt):
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
    A1 = 6.80e8;    E1 = 156      # wood -> gas
    A2 = 8.23e8;    E2 = 148      # wood -> tar
    A3 = 2.91e2;    E3 = 61       # wood -> char
    
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
            
    # return the wood, char, gas, tar concentrations as a density, kg/m^3
    return pw, pg, pt, pc


# Function - primary reactions from pyroprobe 100
# -----------------------------------------------------------------------------

def font2(rhow, T, dt, nt):
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
    A1 = 1.52e7;    E1 = 139      # wood -> gas
    A2 = 5.85e6;    E2 = 119      # wood -> tar
    A3 = 2.98e3;    E3 = 73       # wood -> char
    
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
       
    # return the wood, char, gas, tar concentrations as a density, kg/m^3
    return pw, pg, pt, pc
    