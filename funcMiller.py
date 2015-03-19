"""
Functions based on Miller 1997 kinetic reaction scheme for biomass pyrolysis. 
Reactions are evaluated at some temperature.

Reference:
Miller, Bellan, 1997. Combust. Sci. and Tech., 126, 97-137.
"""

import numpy as np

# Function for cellulose reactions
# -----------------------------------------------------------------------------

def cell(T, pw, dt, p):
    """
    T = temperature, K
    pw = vector of initial wood concentration, kg/m^3
    dt = time step, s
    p = total number of time steps
    """
    
    # vector of initial cellulose concentration for beech wood, Table 2, kg/m^3
    cell = pw*0.48  # cellulose
    
    # vectors to store product concentrations, kg/m^3
    cella = np.zeros(p)     # active cellulose
    gas = np.zeros(p)       # gas
    tar = np.zeros(p)       # tar
    char = np.zeros(p)      # char
    
    R = 0.008314    # universal gas constant, kJ/mol*K
    x = 0.35        # char formation mass ratio for cellulose
    
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 2.8e19;    E1 = 242.4    # cell -> cella
    A2 = 3.28e14;   E2 = 196.5    # cella -> tar
    A3 = 1.3e10;    E3 = 150.5    # cella -> x*char + (1-x)*gas
    A4 = 4.28e6;    E4 = 108      # tar -> gas
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # cell -> cella
    K2 = A2 * np.exp(-E2 / (R * T))  # cella -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # cella -> x*char + (1-x)*gas
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas
    
    # calculate concentrations for each product, kg/m^3
    for i in range(1, p):
        r1 = K1 * cell[i-1]     # cell -> cella
        r2 = K2 * cella[i-1]    # cella -> tar
        r3 = K3 * cella[i-1]    # cella -> x*char + (1-x)*gas
        r4 = K4 * tar[i-1]      # tar -> gas
        cell[i] = cell[i-1] - r1*dt                 # cellulose
        cella[i] = cella[i-1] + r1*dt - (r2+r3)*dt  # active cellulose
        gas[i] = gas[i-1] + r3*(1-x)*dt + r4*dt     # gas
        tar[i] = tar[i-1] + r2*dt - r4*dt           # tar
        char[i] = char[i-1] + r3*x*dt               # char
        
    # return vector for each concentration, kg/m^3
    return cell, cella, gas, tar, char
    
    
# Function for hemicellulose reactions
# -----------------------------------------------------------------------------

def hemi(T, pw, dt, p):
    """
    T = temperature, K
    pw = vector of initial wood concentration, kg/m^3
    dt = time step, s
    p = total number of time steps
    """
    
    # vector of initial hemicellulose concentration for beech wood, see Table 2
    hemi = pw*0.28  # hemicellulose, kg/m^3
    
    # vectors to store product concentrations, kg/m^3
    hemia = np.zeros(p)     # active hemicellulose
    gas = np.zeros(p)       # gas
    tar = np.zeros(p)       # tar
    char = np.zeros(p)      # char
    
    R = 0.008314    # universal gas constant, kJ/mol*K
    x = 0.60        # char formation mass ratio for cellulose
    
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 2.1e16;    E1 = 186.7    # hemi -> hemia
    A2 = 8.75e15;   E2 = 202.4    # hemia -> tar
    A3 = 2.6e11;    E3 = 145.7    # hemia -> x*char + (1-x)*gas
    A4 = 4.28e6;    E4 = 108      # tar -> gas
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # hemi -> hemia
    K2 = A2 * np.exp(-E2 / (R * T))  # hemia -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # hemia -> x*char + (1-x)*gas
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas
    
    # calculate concentrations for each product, kg/m^3
    for i in range(1, p):
        r1 = K1 * hemi[i-1]     # hemi -> hemia
        r2 = K2 * hemia[i-1]    # hemia -> tar
        r3 = K3 * hemia[i-1]    # hemia -> x*char + (1-x)*gas
        r4 = K4 * tar[i-1]      # tar -> gas
        hemi[i] = hemi[i-1] - r1*dt                 # cellulose
        hemia[i] = hemia[i-1] + r1*dt - (r2+r3)*dt  # active cellulose
        gas[i] = gas[i-1] + r3*(1-x)*dt + r4*dt     # gas
        tar[i] = tar[i-1] + r2*dt - r4*dt           # tar
        char[i] = char[i-1] + r3*x*dt               # char
        
    # return vector for each concentration, kg/m^3
    return hemi, hemia, gas, tar, char
    
    
# Function for lignin reactions
# -----------------------------------------------------------------------------

def lig(T, pw, dt, p):
    """
    T = temperature, K
    pw = vector of initial wood concentration, kg/m^3
    dt = time step, s
    p = total number of time steps
    """
    
    # vector of initial lignin concentration for beech wood, see Table 2
    lig = pw*0.24  # hemicellulose, kg/m^3
    
    # vectors to store product concentrations, kg/m^3
    liga = np.zeros(p)      # active lignin
    gas = np.zeros(p)       # gas
    tar = np.zeros(p)       # tar
    char = np.zeros(p)      # char
    
    R = 0.008314    # universal gas constant, kJ/mol*K
    x = 0.75        # char formation mass ratio for cellulose
    
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 2.1e16;    E1 = 186.7    # lig -> liga
    A2 = 8.75e15;   E2 = 202.4    # liga -> tar
    A3 = 2.6e11;    E3 = 145.7    # liga -> x*char + (1-x)*gas
    A4 = 4.28e6;    E4 = 108      # tar -> gas
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # lig -> liga
    K2 = A2 * np.exp(-E2 / (R * T))  # liga -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # liga -> x*char + (1-x)*gas
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas
    
    # calculate concentrations for each product, kg/m^3
    for i in range(1, p):
        r1 = K1 * lig[i-1]     # lig -> liga
        r2 = K2 * liga[i-1]    # liga -> tar
        r3 = K3 * liga[i-1]    # liga -> x*char + (1-x)*gas
        r4 = K4 * tar[i-1]      # tar -> gas
        lig[i] = lig[i-1] - r1*dt                   # lignin
        liga[i] = liga[i-1] + r1*dt - (r2+r3)*dt    # active lignin
        gas[i] = gas[i-1] + r3*(1-x)*dt + r4*dt     # gas
        tar[i] = tar[i-1] + r2*dt - r4*dt           # tar
        char[i] = char[i-1] + r3*x*dt               # char
        
    # return vector for each concentration, kg/m^3
    return lig, liga, gas, tar, char
    
    