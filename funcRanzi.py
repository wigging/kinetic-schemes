"""
Function based on Ranzi 2014 kinetic reaction scheme for biomass pyrolysis. 
Reactions are evaluated at some temperature.

Reference:
Ranzi, Corbetta, Manenti, Pierucci, 2014. Chemical Engineering Science, 110, 2-12.
"""

# Modules
# -----------------------------------------------------------------------------

import numpy as np

# Function - primary kinetic reactions from Table 1 in Supplemental Material
# -----------------------------------------------------------------------------

def ranzi(T, pcell, pcella, plvg, dt, i):
    """
    Primary kinetic reactions
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
    
    R = 1.987   # universal gas constant, kcal/kmol*K
    
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    A1 = 4e13;  E1 = 45000  # CELL -> CELLA
    A2 = 0.5e9; E2 = 29000  # CELLA -> products
    A3 = 1.8;   E3 = 10000  # CELLA -> LVG
    A4 = 4e7;   E4 = 31000  # CELL -> 5*H2O + 6*Char
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # CELL -> CELLA
    K2 = A2 * np.exp(-E2 / (R * T))  # CELLA -> products
    K3 = A3 * T * np.exp(-E3 / (R * T))  # CELLA -> LVG
    K4 = A4 * np.exp(-E4 / (R * T))  # CELL -> 5*H2O + 6*Char
    
    # reaction rate for each reaction, rho/s
    rcell = -(K1+K4) * pcell[i-1]
    rcella = K1 * pcell[i-1] - (K2+K3) * pcella[i-1]
    rlvg = K3 * pcella[i-1]
    
    # species concentrations
    cell = pcell[i-1] + rcell*dt                     # CELL
    cella = pcella[i-1] + rcella*dt
    lvg = plvg[i-1] + rlvg*dt
    
    return cell, cella, lvg
    
    
# Function for cellulose reactions
# -----------------------------------------------------------------------------

def cell(T, pw, dt, p):
    """
    T = temperature, K
    pw = vector of initial wood concentration, kg/m^3
    dt = time step, s
    p = total number of time steps
    """
    
    # array to store species concentrations as a density, kg/m^3
    # row = chemical species
    # column = concentration at time step
    spec = np.zeros([16, p])
    
    # vector of initial cellulose concentration in wood, kg/m^3
    spec[0] = pw*0.5
    
    R = 1.987   # universal gas constant, kcal/kmol*K
    
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    A1 = 4e13;  E1 = 45000  # CELL -> CELLA
    A2 = 0.5e9; E2 = 29000  # CELLA -> products
    A3 = 1.8;   E3 = 10000  # CELLA -> LVG
    A4 = 4e7;   E4 = 31000  # CELL -> 5*H2O + 6*Char
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))         # CELL -> CELLA
    K2 = A2 * np.exp(-E2 / (R * T))         # CELLA -> G2
    K3 = A3 * T * np.exp(-E3 / (R * T))     # CELLA -> LVG
    K4 = A4 * np.exp(-E4 / (R * T))         # CELL -> G1
    
    # calculate concentrations for each chemical species. kg/m^3
#    for i in range(1, p):
#        spec[0, i] = spec[0, i-1] - (K1+K4)*spec[0, i-1]*dt                        # CELL
#        spec[1, i] = spec[1, i-1] + K1*spec[0, i-1]*dt - (K2+K3)*spec[1, i-1]*dt   # CELLA
#        spec[2, i] = spec[2, i-1] + K3*spec[1, i-1]*dt                             # LVG
#        spec[3, i] = spec[3, i-1] + K4*spec[0, i-1]*dt*5 + K2*spec[1, i-1]*dt*0.83 # H2O
#        spec[4, i] = spec[4, i-1] + K2*spec[1, i-1]*dt*0.8                         # HAA
#        spec[5, i] = spec[5, i-1] + K2*spec[1, i-1]*dt*0.2                         # GLYOX
#        spec[6, i] = spec[6, i-1] + K2*spec[1, i-1]*dt*0.1                         # C2H4O
#        spec[7, i] = spec[7, i-1] + K2*spec[1, i-1]*dt*0.25                        # HMFU
#        spec[8, i] = spec[8, i-1] + K2*spec[1, i-1]*dt*0.3                         # C3H6O
#        spec[9, i] = spec[9, i-1] + K2*spec[1, i-1]*dt*0.21                        # CO2
#        spec[10, i] = spec[10, i-1] + K2*spec[1, i-1]*dt*0.1                       # H2
#        spec[11, i] = spec[11, i-1] + K2*spec[1, i-1]*dt*0.4                       # CH2O
#        spec[12, i] = spec[12, i-1] + K2*spec[1, i-1]*dt*0.16                      # CO
#        spec[13, i] = spec[13, i-1] + K2*spec[1, i-1]*dt*0.1                       # CH4
#        spec[14, i] = spec[14, i-1] + K2*spec[1, i-1]*dt*0.02                      # HCOOH
#        spec[15, i] = spec[15, i-1] + K4*spec[0, i-1]*dt*6 + K2*spec[1, i-1]*dt*0.61 # Char
    
    # concentrations
    for i in range(1, p):
        r1 = K1 * spec[0, i-1]
        r2 = K2 * spec[1, i-1]
        r3 = K3 * spec[1, i-1]
        r4 = K4 * spec[0, i-1]
        spec[0, i] = spec[0, i-1] - (r1+r4)*dt          # CELL
        spec[1, i] = spec[1, i-1] + r1*dt - (r2+r3)*dt  # CELLA
        spec[2, i] = spec[2, i-1] + r4*dt               # G1
        spec[3, i] = spec[3, i-1] + r2*dt               # G2
        spec[4, i] = spec[4, i-1] + r3*dt               # LVG
    
    # return species array concentrations, kg/m^3
    return spec
    
    
# Function for cellulose reactions
# -----------------------------------------------------------------------------

def cell2(T, pw, dt, p):
    """
    T = temperature, K
    pw = vector of initial wood concentration, kg/m^3
    dt = time step, s
    p = total number of time steps
    """
    
    # vector of initial cellulose concentration in wood, kg/m^3
    cell = pw*0.5
    
    # vectors to store product concentrations, kg/m^3
    cella = np.zeros(p)
    g1 = np.zeros(p)
    g2 = np.zeros(p)
    lvg = np.zeros(p)
    
    R = 1.987   # universal gas constant, kcal/kmol*K
    
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    A1 = 4e13;  E1 = 45000  # CELL -> CELLA
    A2 = 0.5e9; E2 = 29000  # CELLA -> G2
    A3 = 1.8;   E3 = 10000  # CELLA -> LVG
    A4 = 4e7;   E4 = 31000  # CELL -> G1
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))         # CELL -> CELLA
    K2 = A2 * np.exp(-E2 / (R * T))         # CELLA -> G2
    K3 = A3 * T * np.exp(-E3 / (R * T))     # CELLA -> LVG
    K4 = A4 * np.exp(-E4 / (R * T))         # CELL -> G1
    
    # calculate concentrations for each product, kg/m^3
    
    # concentrations
    for i in range(1, p):
        r1 = K1 * cell[i-1]
        r2 = K2 * cella[i-1]
        r3 = K3 * cella[i-1]
        r4 = K4 * cell[i-1]
        cell[i] = cell[i-1] - (r1+r4)*dt            # CELL
        cella[i] = cella[i-1] + r1*dt - (r2+r3)*dt  # CELLA
        g1[i] = g1[i-1] + r4*dt                     # G1
        g2[i] = g2[i-1] + r2*dt                     # G2
        lvg[i] = lvg[i-1] + r3*dt                   # LVG
    
    # return species array concentrations, kg/m^3
    return cell, cella, g1, g2, lvg


# Function for cellulose reactions
# -----------------------------------------------------------------------------

def cell3(T, mw, dt, p):
    """
    T = temperature, K
    pw = vector of initial wood concentration, kg/m^3
    dt = time step, s
    p = total number of time steps
    """
    
    # vector of initial cellulose concentration in wood, kg/m^3
    # cell = pw*0.5
    