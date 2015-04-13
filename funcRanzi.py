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
    
    # vectors to store main product concentrations, kg/m^3
    cell = pw*(wt/100)      # initial cellulose conc. in wood
    g1 = np.zeros(nt)       # G1
    cella = np.zeros(nt)    # CELLA
    lvg = np.zeros(nt)      # LVG
    g4 = np.zeros(nt)       # G4
    
    R = 1.987   # universal gas constant, kcal/kmol*K
    
    # reaction rate constant for each reaction, 1/s
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    K1 = 4e7 * np.exp(-31000 / (R * T))         # CELL -> G1
    K2 = 4e13 * np.exp(-45000 / (R * T))        # CELL -> CELLA
    K3 = 1.8 * T * np.exp(-10000 / (R * T))     # CELLA -> LVG
    K4 = 0.5e9 * np.exp(-29000 / (R * T))       # CELLA -> G4
    
    # sum of moles in each group, mol
    sumg1 = 11      # sum of G1
    sumg4 = 4.08    # sum of G4
    
    # calculate concentrations for main groups, kg/m^3
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
    
    # store main groups in array
    main = np.array([cell, g1, cella, lvg, g4])
    
    # total group concentration per total moles in that group, (kg/m^3) / mol
    fg1 = g1/sumg1  # fraction of G1 
    fg4 = g4/sumg4  # fraction of G4
    
    # array to store product concentrations as a density, kg/m^3
    prod = np.zeros([21, nt])
    prod[0] = 0.16*fg4              # CO
    prod[1] = 0.21*fg4              # CO2
    prod[2] = 0.4*fg4               # CH2O
    prod[3] = 0.02*fg4              # HCOOH
    prod[5] = 0.1*fg4               # CH4
    prod[6] = 0.2*fg4               # Glyox
    prod[8] = 0.1*fg4               # C2H4O
    prod[9] = 0.8*fg4               # HAA
    prod[11] = 0.3*fg4              # C3H6O
    prod[14] = 0.25*fg4             # HMFU
    prod[15] = lvg                  # LVG
    prod[18] = 0.1*fg4              # H2
    prod[19] = 5*fg1 + 0.83*fg4     # H2O
    prod[20] = 6*fg1 + 0.61*fg4     # Char
    
    # return arrays of main groups and products as mass fraction, (-)
    return main/rhow, prod/rhow
    
    
# Function - hemicellulose reactions, HCE
# -----------------------------------------------------------------------------

def hemi(rhow, wt, T, dt, nt):
    """
    rhow = wood density, kg/m^3
    wt = weight percent of hemicellulose, %
    T = temperature, K
    dt = time step, s
    nt = total number of time steps
    """
    
    # vector for initial wood concentration, kg/m^3
    pw = np.ones(nt)*rhow
    
    # vectors to store main product concentrations, kg/m^3
    hce = pw*(wt/100)   # initial hemicellulose conc. in wood
    g1 = np.zeros(nt)   # G1
    g2 = np.zeros(nt)   # G2
    g3 = np.zeros(nt)   # G3
    g4 = np.zeros(nt)   # G4
    xyl = np.zeros(nt)  # Xylan
    
    R = 1.987   # universal gas constant, kcal/kmol*K
    
    # reaction rate constant for each reaction, 1/s
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    K1 = 0.33e10 * np.exp(-31000 / (R * T))     # HCE -> G1
    K2 = 0.33e10 * np.exp(-33000 / (R * T))     # HCE2 -> G2
    K3 = 0.05 * T * np.exp(-8000 / (R * T))     # HCE1 -> G3
    K4 = 1e9 * np.exp(-32000 / (R * T))         # HCE1 -> G4
    K5 = 0.9 * T * np.exp(-11000 / (R * T))     # HCE1 -> Xylan
    
    # sum of moles in each group, mol
    sumg2 = 4.625   # sum of G2
    sumg3 = 4.875   # sum of G3
    sumg4 = 4.775   # sum of G4
    
    # calculate concentrations for main groups, kg/m^3
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
        
    # store main groups in array
    main = np.array([hce, g1, g2, g3, g4, xyl])
    
    # total group concentration per total moles in that group, (kg/m^3)/mol
    fg2 = g2/sumg2  # fraction of G2 
    fg3 = g3/sumg3  # fraction of G3
    fg4 = g4/sumg4  # fraction of G4
    
    # array to store product concentrations as a density, kg/m^3
    prod = np.zeros([21, nt])
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
    
    # return arrays of main groups and products as mass fraction, (-)
    return main/rhow, prod/rhow
    
    
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
    
    # vectors to store main product concentrations, kg/m^3
    ligc = pw*(wt/100/3)    # initial lignin in wood, assume 1/3 of total lignin
    g1 = np.zeros(nt)
    g2 = np.zeros(nt)
    
    R = 1.987       # universal gas constant, kcal/kmol*K
    
    # reaction rate constant for each reaction, 1/s
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    K1 = 1.33e15 * np.exp(-48500 / (R * T))     # LIG-C -> G1
    K2 = 1.6e6 * np.exp(-31500 / (R * T))       # LIG-CC -> G2
    
    # sum of moles in each group, mol
    sumg1 = 9.49    # sum of G1
    sumg2 = 11.35   # sum of G2

    # calculate concentrations for main groups, kg/m^3
    for i in range(1, nt):
        r1 = K1 * ligc[i-1]             # LIG-C -> G1
        r2 = K2 * 0.35*g1[i-1]/sumg1    # LIG-CC -> G2
        ligc[i] = ligc[i-1] - r1*dt         # LIG-C
        g1[i] = g1[i-1] + r1*dt - r2*dt     # G1
        g2[i] = g2[i-1] + r2*dt             # G2
    
    # store main groups in array
    main = np.array([ligc, g1, g2])
    
    # total group concentration per total moles in that group, (kg/m^3)/mol
    fg1 = g1/sumg1  # fraction of G1 
    fg2 = g2/sumg2  # fraction of G2
    
    # array to store product concentrations as a density, kg/m^3
    prod = np.zeros([21, nt])
    prod[0] = 0.32*fg1 + (0.4 + 0.4)*fg2    # CO
    prod[2] = (0.3 + 0.7)*fg1 + 1*fg2       # CH2O
    prod[5] = 0.495*fg1 + 0.65*fg2          # CH4
    prod[7] = 0.41*fg1 + 0.6*fg2            # C2H4
    prod[9] = 0.35*fg2                      # HAA
    prod[13] = 0.08*fg1 + 0.2*fg2           # Phenol
    prod[16] = 0.1*fg1 + 0.3*fg2            # Coumaryl
    prod[19] = 1*fg1 + 0.7*fg2              # H2O
    prod[20] = 5.735*fg1 + 6.75*fg2         # Char
    
    # return arrays of main groups and products as mass fractions, (-)
    return main/rhow, prod/rhow
    
    
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
    
    # vectors to store main product concentrations, kg/m^3
    ligh = pw*(wt/100/3)    # initial lignin in wood, assume 1/3 of total lignin
    g1 = np.zeros(nt)       # G1
    g2 = np.zeros(nt)       # G2
    g3 = np.zeros(nt)       # G3
    g4 = np.zeros(nt)       # G4
    g5 = np.zeros(nt)       # G4
    fe2macr = np.zeros(nt)  # FE2MACR
    
    R = 1.987       # universal gas constant, kcal/kmol*K   
    
    # reaction rate constant for each reaction, 1/s
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    K1 = 0.67e13 * np.exp(-37500 / (R * T))     # LIG-H -> G1
    K2 = 33 * np.exp(-15000 / (R * T))          # LIG-OH -> G2
    K3 = 0.5e8 * np.exp(-30000 / (R * T))       # LIG-OH -> LIG
    K4 = 0.083 * T * np.exp(-8000 / (R * T))    # LIG -> G4
    K5 = 0.4e9 * np.exp(-30000 / (R * T))       # LIG -> G5
    K6 = 2.4 * T * np.exp(-12000 / (R * T))     # LIG -> FE2MACR
    
    # sum of moles in each group, mol
    sumg1 = 2       # sum of G1
    sumg2 = 20.7    # sum of G2
    sumg3 = 9.85    # sum of G3
    sumg4 = 11.1    # sum of G4
    sumg5 = 10.7    # sum of G5
    
    # calculate concentrations for main groups, kg/m^3
    for i in range(1, nt):
        r1 = K1 * ligh[i-1]             # LIG-H -> G1
        r2 = K2 * 1*g1[i-1]/sumg1       # LIG-OH -> G2
        r3 = K3 * 1*g1[i-1]/sumg1       # LIG-OH -> LIG
        r4 = K4 * 1*g3[i-1]/sumg3       # LIG -> G4
        r5 = K5 * 1*g3[i-1]/sumg3       # LIG -> G5
        r6 = K6 * 1*g3[i-1]/sumg3       # LIG -> FE2MACR
        ligh[i] = ligh[i-1] - r1*dt                 # LIG-H
        g1[i] = g1[i-1] + r1*dt - (r2+r3)*dt        # G1
        g2[i] = g2[i-1] + r2*dt                     # G2
        g3[i] = g3[i-1] + r3*dt - (r4+r5+r6)*dt     # G3
        g4[i] = g4[i-1] + r4*dt                     # G4
        g5[i] = g5[i-1] + r5*dt                     # G5
        fe2macr[i] = fe2macr[i-1] + r6*dt           # FE2MACR
    
    # store main groups in array
    main = np.array([ligh, g1, g2, g3, g4, g5, fe2macr])
    
    # total group concentration per total moles in that group, (kg/m^3)/mol
    fg1 = g1/sumg1  # fraction of G1 
    fg2 = g2/sumg2  # fraction of G2
    fg3 = g3/sumg3  # fraction of G3
    fg4 = g4/sumg4  # fraction of G4
    fg5 = g5/sumg5  # fraction of G5
    
    # array to store product concentrations as a density, kg/m^3
    prod = np.zeros([21, nt])
    prod[0] = (0.5 + 1.6)*fg2 +(0.3 + 1)*fg3 +(0.4 + 0.2)*fg4 +(1 + 0.45)*fg5 # CO
    prod[1] = 0.05*fg3                                              # CO2
    prod[2] = 3.9*fg2 + 0.6*fg3 + (2 + 0.4)*fg4 + (0.2 + 0.5)*fg5   # CH2O
    prod[3] = 0.05*fg3 + 0.05*fg5                                   # HCOOH
    prod[4] = 0.5*fg2 + (0.5 + 0.5)*fg3 + 0.4*fg4 + 0.4*fg5         # CH3OH
    prod[5] = (0.1 + 1.65)*fg2 +(0.1 + 0.35)*fg3 +(0.2 + 0.4)*fg4 +(0.2 + 0.4)*fg5 # CH4
    prod[6] = 0                                                     # Glyox
    prod[7] = 0.3*fg2 + 0.2*fg3 + 0.5*fg4 + 0.65*fg5                # C2H4
    prod[8] = 0.2*fg5                                               # C2H4O
    prod[9] = 0                                                     # HAA
    prod[10] = 0                                                    # C2H5OH
    prod[11] = 1*fg1 + 0.2*fg5                                      # C3H6O
    prod[12] = 0                                                    # Xylan
    prod[13] = 0                                                    # Phenol
    prod[14] = 0                                                    # HMFU
    prod[15] = 0                                                    # LVG
    prod[16] = 0                                                    # Coumaryl
    prod[17] = fe2macr                                              # FE2MACR
    prod[18] = 0.5*fg2 + 0.15*fg3                                   # H2
    prod[19] = 1.5*fg2 + 0.9*fg3 + 0.6*fg4 + 0.95*fg5               # H2O
    prod[20] = 10.15*fg2 + 4.15*fg3 + 6*fg4 + 5.5*fg5               # Char
    
    # return arrays of main groups and products as mass fractions, (-)
    return main/rhow, prod/rhow
    
    
# Function - lignin oxygen-rich reactions, LIG-O
# -----------------------------------------------------------------------------

def ligo(rhow, wt, T, dt, nt):
    """
    rhow = wood density, kg/m^3
    wt = weight percent of lignin-h, %
    T = temperature, K
    dt = time step, s
    nt = total number of time steps
    """
    
    # vector for initial wood concentration, kg/m^3
    pw = np.ones(nt)*rhow
    
    # vectors to store main product concentrations, kg/m^3
    ligo = pw*(wt/100/3)    # initial lignin in wood, assume 1/3 of total lignin
    g1 = np.zeros(nt)       # G1
    g2 = np.zeros(nt)       # G2
    g3 = np.zeros(nt)       # G3
    g4 = np.zeros(nt)       # G4
    g5 = np.zeros(nt)       # G4
    fe2macr = np.zeros(nt)  # FE2MACR
    
    R = 1.987       # universal gas constant, kcal/kmol*K
    
    # reaction rate constant for each reaction, 1/s
    # A = pre-factor (1/s) and E = activation energy (kcal/kmol)
    K1 = 0.33e9 * np.exp(-25500 / (R * T))      # LIG-O -> G1
    K2 = 33 * np.exp(-15000 / (R * T))          # LIG-OH -> G2
    K3 = 0.5e8 * np.exp(-30000 / (R * T))       # LIG-OH -> LIG
    K4 = 0.083 * T * np.exp(-8000 / (R * T))    # LIG -> G4
    K5 = 0.4e9 * np.exp(-30000 / (R * T))       # LIG -> G5
    K6 = 2.4 * T * np.exp(-12000 / (R * T))     # LIG -> FE2MACR
    
    # sum of moles in each group, mol
    sumg1 = 2       # sum of G1
    sumg2 = 20.7    # sum of G2
    sumg3 = 9.85    # sum of G3
    sumg4 = 11.1    # sum of G4
    sumg5 = 10.7    # sum of G5
    
    # calculate concentrations for main groups, kg/m^3
    for i in range(1, nt):
        r1 = K1 * ligo[i-1]             # LIG-O -> G1
        r2 = K2 * 1*g1[i-1]/sumg1       # LIG-OH -> G2
        r3 = K3 * 1*g1[i-1]/sumg1       # LIG-OH -> LIG
        r4 = K4 * 1*g3[i-1]/sumg3       # LIG -> G4
        r5 = K5 * 1*g3[i-1]/sumg3       # LIG -> G5
        r6 = K6 * 1*g3[i-1]/sumg3       # LIG -> FE2MACR
        ligo[i] = ligo[i-1] - r1*dt                 # LIG-H
        g1[i] = g1[i-1] + r1*dt - (r2+r3)*dt        # G1
        g2[i] = g2[i-1] + r2*dt                     # G2
        g3[i] = g3[i-1] + r3*dt - (r4+r5+r6)*dt     # G3
        g4[i] = g4[i-1] + r4*dt                     # G4
        g5[i] = g5[i-1] + r5*dt                     # G5
        fe2macr[i] = fe2macr[i-1] + r6*dt           # FE2MACR
    
    # store main groups in array
    main = np.array([ligo, g1, g2, g3, g4, g5, fe2macr])
    
    # total group concentration per total moles in that group, (kg/m^3)/mol
    fg1 = g1/sumg1  # fraction of G1 
    fg2 = g2/sumg2  # fraction of G2
    fg3 = g3/sumg3  # fraction of G3
    fg4 = g4/sumg4  # fraction of G4
    fg5 = g5/sumg5  # fraction of G5
    
    # array to store product concentrations as a density, kg/m^3
    prod = np.zeros([21, nt])
    prod[0] = (0.5 + 1.6)*fg2 +(0.3 + 1)*fg3 +(0.4 + 0.2)*fg4 +(1 + 0.45)*fg5 # CO
    prod[1] = 1*fg1 + 0.05*fg3                                      # CO2
    prod[2] = 3.9*fg2 + 0.6*fg3 + (2 + 0.4)*fg4 + (0.2 + 0.5)*fg5   # CH2O
    prod[3] = 0.05*fg3 + 0.05*fg5                                   # HCOOH
    prod[4] = 0.5*fg2 + (0.5 + 0.5)*fg3 + 0.4*fg4 + 0.4*fg5         # CH3OH
    prod[5] = (0.1 + 1.65)*fg2 +(0.1 + 0.35)*fg3 +(0.2 + 0.4)*fg4 +(0.2 + 0.4)*fg5 # CH4
    prod[6] = 0                                                     # Glyox
    prod[7] = 0.3*fg2 + 0.2*fg3 + 0.5*fg4 + 0.65*fg5                # C2H4
    prod[8] = 0.2*fg5                                               # C2H4O
    prod[9] = 0                                                     # HAA
    prod[10] = 0                                                    # C2H5OH
    prod[11] = 0.2*fg5                                              # C3H6O
    prod[12] = 0                                                    # Xylan
    prod[13] = 0                                                    # Phenol
    prod[14] = 0                                                    # HMFU
    prod[15] = 0                                                    # LVG
    prod[16] = 0                                                    # Coumaryl
    prod[17] = fe2macr                                              # FE2MACR
    prod[18] = 0.5*fg2 + 0.15*fg3                                   # H2
    prod[19] = 1.5*fg2 + 0.9*fg3 + 0.6*fg4 + 0.95*fg5               # H2O
    prod[20] = 10.15*fg2 + 4.15*fg3 + 6*fg4 + 5.5*fg5               # Char
    
    # return arrays of main groups and products as mass fractions, (-)
    return main/rhow, prod/rhow
    
    