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
    
    # vectors to store product concentrations, kg/m^3
    cell = pw*(wt/100)      # initial cellulose conc. in wood
    g1 = np.zeros(nt)       # G1
    cella = np.zeros(nt)    # CELLA
    lvg = np.zeros(nt)      # LVG
    g4 = np.zeros(nt)       # G4
    
    R = 1.987   # universal gas constant, kcal/kmol*K
    
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    A1 = 4e7;   E1 = 31000      # CELL -> G1
    A2 = 4e13;  E2 = 45000      # CELL -> CELLA
    A3 = 1.8*T;   E3 = 10000    # CELLA -> LVG
    A4 = 0.5e9; E4 = 29000      # CELLA -> G4
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))     # CELL -> G1
    K2 = A2 * np.exp(-E2 / (R * T))     # CELL -> CELLA
    K3 = A3 * np.exp(-E3 / (R * T))     # CELLA -> LVG
    K4 = A4 * np.exp(-E4 / (R * T))     # CELLA -> G4
    
    # calculate concentrations for each product, kg/m^3
    for i in range(1, nt):
        r1 = K1 * cell[i-1]     # CELL -> G1
        r2 = K2 * cell[i-1]     # CELL -> CELLA
        r3 = K3 * cella[i-1]    # CELLA -> LVG
        r4 = K4 * cella[i-1]    # CELLA -> G4
        cell[i] = cell[i-1] - (r1+r2)*dt            # CELL
        g1[i] = g1[i-1] + r1*dt                     # G1
        cella[i] = cella[i-1] + r2*dt - (r3+r4)*dt  # CELLA
        lvg[i] = lvg[i-1] + r3*dt                   # LVG
        g4[i] = g4[i-1] + r4*dt                     # G4
    
    sumg1 = 11      # sum of G1 species
    sumg4 = 4.08    # sum of G4 species
    
    h2o = 5*g1/sumg1 + 0.83*g4/sumg4    # H2O
    char = 6*g1/sumg1 + 0.61*g4/sumg4   # Char
    haa = 0.8*g4/sumg4                  # HAA
    glyox = 0.2*g4/sumg4                # Glyox
    c2h4o = 0.1*g4/sumg4                # C2H4O
    hmfu = 0.25*g4/sumg4                # HMFU
    c3h6o = 0.3*g4/sumg4                # C3H6O
    co2 = 0.21*g4/sumg4                 # CO2
    h2 = 0.1*g4/sumg4                   # H2
    ch2o = 0.4*g4/sumg4                 # CH2O
    co = 0.16*g4/sumg4                  # CO
    ch4 = 0.1*g4/sumg4                  # CH4
    hcooh = 0.02*g4/sumg4               # HCOOH
    
    # return array of concentrations as a density, kg/m^3
    prod = np.array([cell, cella, g1, g4, lvg, h2o, char, haa, glyox, c2h4o, \
                    hmfu, c3h6o, co2, h2, ch2o, co, ch4, hcooh])
    
    return prod
    
    
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
    
    # vectors to store product concentrations, kg/m^3
    hce = pw*(wt/100)   # initial hemicellulose conc. in wood
    g1 = np.zeros(nt)   # G1
    g2 = np.zeros(nt)   # G2
    g3 = np.zeros(nt)   # G3
    g4 = np.zeros(nt)   # G4
    xyl = np.zeros(nt)  # Xylan
    
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
    # where HCE1 as 0.4*g1/(0.4+0.6) and HCE2 as 0.6*g1/(0.4+0.6)
    for i in range(1, nt):
        r1 = K1 * hce[i-1]      # HCE -> G1
        r2 = K2 * 0.6*g1[i-1]   # HCE2 -> G2
        r3 = K3 * 0.4*g1[i-1]   # HCE1 -> G3
        r4 = K4 * 0.4*g1[i-1]   # HCE1 -> G4
        r5 = K5 * 0.4*g1[i-1]   # HCE1 -> Xylan
        hce[i] = hce[i-1] - r1*dt                   # HCE
        g1[i] = g1[i-1] + r1*dt - (r2+r3+r4+r5)*dt  # G1
        g2[i] = g2[i-1] + r2*dt                     # G2
        g3[i] = g3[i-1] + r3*dt                     # G3
        g4[i] = g4[i-1] + r4*dt                     # G4
        xyl[i] = xyl[i-1] + r5*dt                   # Xylan
        
    # sum of moles in each group, mol
    sumg2 = 4.625   # sum of G2
    sumg3 = 4.875   # sum of G3
    sumg4 = 4.775   # sum of G4
    
    # group concentration per total moles in that group, (kg/m^3)/mol
    fg2 = g2/sumg2  # fraction of G2 
    fg3 = g3/sumg3  # fraction of G3
    fg4 = g4/sumg4  # fraction of G4
    
    prod = np.zeros([21, nt]) # array to store products
    
    prod[0] = 0.175*fg2 + (0.3 + 0.15)*fg3 + 0.5*fg4                # CO
    prod[1] = (0.275+0.4)*fg2 + (0.5+0.25)*fg3 + (0.5+0.275)*fg4    # CO2
    prod[2] = (0.5+0.925)*fg2 + 1.7*fg3 + (0.8+0.4)*fg4             # CH2O
    prod[3] = 0.025*fg2 + 0.05*fg3 + 0.025*fg4                      # HCOOH
    prod[4] = 0.3*fg2 + (0.1+0.45)*fg4                              # CH3OH
    prod[5] = 0.25*fg2 + 0.625*fg3 + 0.325*fg4                      # CH4
    prod[7] = 0.275*fg2 + 0.375*fg3 + 0.25*fg4                      # C2H4
    prod[9] = 0.2*fg2                                               # HAA
    prod[10] = 0.1*fg2 + 0.125*fg4                                  # C2H5OH
    prod[12] = xyl                                                  # Xylan
    prod[18] = 0.125*fg4                                            # H2
    prod[19] = 0.2*fg2 + 0.25*fg3 + 0.025*fg4                       # H2O
    prod[20] = 1*fg2 + 0.675*fg3 + 0.875*fg4                        # Char
    
    # return array of concentrations as a density, kg/m^3
    prod = np.array([hce, g1, g2, g3, g4, xyl])
    
    return prod
    
    
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
    
