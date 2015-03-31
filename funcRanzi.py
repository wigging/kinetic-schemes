"""
Function based on Ranzi 2014 kinetic reaction scheme for biomass pyrolysis. 
Reactions are evaluated at some temperature.

Reference:
Ranzi, Corbetta, Manenti, Pierucci, 2014. Chemical Engineering Science, 110, 2-12.
"""

# Modules
# -----------------------------------------------------------------------------

import numpy as np

# Function - calculate molar concentration of biomass
# -----------------------------------------------------------------------------

def conc(cwt, hwt, owt, nwt, ma, rhow):
    """
    cwt = carbon weight percent, %
    hwt = hydrogen weight percent, %
    owt = oxygen weight percent, %
    nwt = nitrogen weight percent, %
    ma = mass of sample (typically assume 100 grams), g
    rhow = initial concentration of wood as density, kg/m^3
    """
    
    # molar mass of each element, g/mol
    cmm = 12.01; hmm = 1.01; omm = 16; nmm = 14.01
    
    # convert wt % to mass based on sample
    c = ma*(cwt/100)
    h = ma*(hwt/100)
    o = ma*(owt/100)
    n = ma*(nwt/100)
    print('c = {}, h = {}, o = {}, n = {}'.format(c, h, o, n))
    
    # calculate moles of each element in the sample, mol
    cmol = c/cmm
    hmol = h/hmm
    omol = o/omm
    nmol = n/nmm
    print('cmol = {}, hmol = {}, omol = {}, nmol = {}'.format(cmol, hmol, omol, nmol))

    # molecular weight of the sample, g/mol
    mw = cmol*cmm + hmol*hmm + omol*omm + nmol*nmm
    print('mw = ', mw)
    
    # molar concentration, mol/m^3
    conc = (rhow*1000)/mw
    
    # return molar concentration, mol/m^3
    return conc
    
    
# Function - cellulose reactions, CELL
# -----------------------------------------------------------------------------

def cell(rhow, wt, T, dt, nt):
    """
    rhow = wood density, kg/m^3
    wt = weight percent of cellulose, %
    T = temperature, K
    dt = time step, s
    nt = total number of time steps
    """
    
    # vector for initial wood concentration, kg/m^3
    pw = np.ones(nt)*rhow
    
    # vector of initial cellulose concentration in wood, kg/m^3
    cell = pw*(wt/100)
    
    # vectors to store product concentrations, kg/m^3
    cella = np.zeros(nt)
    g1 = np.zeros(nt)
    g2 = np.zeros(nt)
    lvg = np.zeros(nt)
    
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
    for i in range(1, nt):
        r1 = K1 * cell[i-1]
        r2 = K2 * cella[i-1]
        r3 = K3 * cella[i-1]
        r4 = K4 * cell[i-1]
        cell[i] = cell[i-1] - (r1+r4)*dt            # CELL
        cella[i] = cella[i-1] + r1*dt - (r2+r3)*dt  # CELLA
        g1[i] = g1[i-1] + r4*dt                     # G1
        g2[i] = g2[i-1] + r2*dt                     # G2
        lvg[i] = lvg[i-1] + r3*dt                   # LVG
    
    # return array of concentrations as a density, kg/m^3
    return np.array([cell, cella, g1, g2, lvg])
    
    
# Function - hemicellulose reactions, HCE
# -----------------------------------------------------------------------------

def hce(rhow, wt, T, dt, nt):
    """
    rhow = wood density, kg/m^3
    wt = weight percent of hemicellulose, %
    T = temperature, K
    dt = time step, s
    nt = total number of time steps
    """
    
    # vector for initial wood concentration, kg/m^3
    pw = np.ones(nt)*rhow
    
    # vector of initial hemicellulose concentration in wood, kg/m^3
    hce = pw*(wt/100)
    
    # vectors to store product concentrations, kg/m^3
    g1 = np.zeros(nt)
    g2 = np.zeros(nt)
    g3 = np.zeros(nt)
    g4 = np.zeros(nt)
    xyl = np.zeros(nt)
    
    R = 1.987   # universal gas constant, kcal/kmol*K
    
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    A1 = 0.33e10;  E1 = 31000   # HCE -> G1
    A2 = 0.33e10;  E2 = 33000   # HCE2 -> G2
    A3 = 0.05*T;   E3 = 8000    # HCE1 -> G3
    A4 = 1e9;      E4 = 32000   # HCE1 -> G4
    A5 = 0.9*T;    E5 = 11000   # HCE1 -> Xylan
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))     # HCE -> G1
    K2 = A2 * np.exp(-E2 / (R * T))     # HCE2 -> G2
    K3 = A3 * np.exp(-E3 / (R * T))     # HCE1 -> G3
    K4 = A4 * np.exp(-E4 / (R * T))     # HCE1 -> G4
    K5 = A5 * np.exp(-E5 / (R * T))     # HCE1 -> Xylan
    
    # calculate concentrations for each product, kg/m^3
    for i in range(1, nt):
        r1 = K1 * hce[i-1]
        r2 = K2 * 0.6*g1[i-1]   # from HCE2 as 0.6*g1/(0.4+0.6)
        r3 = K3 * 0.4*g1[i-1]   # from HCE1 as 0.4*g1/(0.4+0.6)
        r4 = K4 * 0.4*g1[i-1]   # also from HCE1
        r5 = K5 * 0.4*g1[i-1]   # also from HCE1
        hce[i] = hce[i-1] - r1*dt                 # HCE
        g1[i] = g1[i-1] + r1*dt - (r2+r3+r4+r5)*dt  # G1
        g2[i] = g2[i-1] + r2*dt                     # G2
        g3[i] = g3[i-1] + r3*dt                     # G3
        g4[i] = g4[i-1] + r4*dt                     # G4
        xyl[i] = xyl[i-1] + r5*dt                   # Xylan
        
    # return array of concentrations as a density, kg/m^3
    return np.array([hce, g1, g2, g3, g4, xyl])
    
    
# Function - lignin carbon-rich reactions, LIG-C
# -----------------------------------------------------------------------------

def ligc(rhow, wt, T, dt, nt):
    """
    rhow = wood density, kg/m^3
    wt = weight percent of lignin-c, %
    T = temperature, K
    dt = time step, s
    nt = total number of time steps
    """
    
    # vector for initial wood concentration, kg/m^3
    pw = np.ones(nt)*rhow
    
    # vector of initial lignin concentration in wood, kg/m^3
    ligc = pw*(wt/100/3) # assume 1/3 of total lignin
    
    # vectors to store product concentrations, kg/m^3
    g1 = np.zeros(nt)
    g2 = np.zeros(nt)
    
    R = 1.987       # universal gas constant, kcal/kmol*K
    sumg1 = 9.49    # sum of G1 species

    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    A1 = 1.33e15;   E1 = 48500   # LIG-C -> G1
    A2 = 1.6e6;     E2 = 31500   # LIG-CC -> G2
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))     # LIG-C -> G1
    K2 = A2 * np.exp(-E2 / (R * T))     # LIG-CC -> G2
    
    # calculate concentrations for each product, kg/m^3
    for i in range(1, nt):
        r1 = K1 * ligc[i-1]
        r2 = K2 * 0.35*g1[i-1]/sumg1 # from LIG-CC as 0.35*g1/(sum G1 species)
        ligc[i] = ligc[i-1] - r1*dt
        g1[i] = g1[i-1] + r1*dt - r2*dt
        g2[i] = g2[i-1] + r2*dt
        
    # return array of concentrations as a density, kg/m^3
    return np.array([ligc, g1, g2])
    
    
# Function - lignin hydrogen-rich reactions, LIG-H
# -----------------------------------------------------------------------------

def ligh(rhow, wt, T, dt, nt):
    """
    rhow = wood density, kg/m^3
    wt = weight percent of lignin-h, %
    T = temperature, K
    dt = time step, s
    nt = total number of time steps
    """
    
    # vector for initial wood concentration, kg/m^3
    pw = np.ones(nt)*rhow
    
    # vector of initial lignin concentration in wood, kg/m^3
    ligh = pw*0.2    
    
