"""
Plot gas, tar, char from primary and secondary reactions as determined by the
Miller and Bellan 1997 kinetic scheme for biomass pyrolysis. Biomass composition
by mass of beech wood is used from Table 2 in paper. This version ignores the
first reaction pathway from virgin -> active to allow larger time steps.

Reference:
Miller, Bellan, 1997. Combust. Sci. and Tech., 126, pp 97-137.
"""

import numpy as np
import matplotlib.pyplot as py

# Parameters
# ------------------------------------------------------------------------------

T = 773  # ambient temp, K

# composition of beech wood from Table 2 in paper
wtcell = 0.48   # cellulose mass fraction, (-)
wthemi = 0.28   # hemicellulose mass fraction, (-)
wtlig = 0.24    # lignin mass fraction, (-)

dt = 0.005                              # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Functions for Miller and Bellan 1997 Kinetic Scheme
# ------------------------------------------------------------------------------

def millercell_noR1(cella, gas, tar, char, T, dt, s=1):
    """
    Cellulose reactions from Miller and Bellan 1997 paper for biomass pyrolysis.

    Parameters
    ----------
    cella = activie cellulose concentration, kg/m^3
    gas = gas concentration, kg/m^3
    tar = tar concentration, kg/m^3
    char = char concentration, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    ncella = new active cellulose concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A2 = 3.28e14;   E2 = 196.5    # cella -> tar
    A3 = 1.3e10;    E3 = 150.5    # cella -> x*char + (1-x)*gas
    A4 = 4.28e6;    E4 = 108      # tar -> gas
    R = 0.008314    # universal gas constant, kJ/mol*K
    x = 0.35        # char formation mass ratio for cellulose

    # reaction rate constant for each reaction, 1/s
    K2 = A2 * np.exp(-E2 / (R * T))  # cella -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # cella -> x*char + (1-x)*gas
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas

    if s == 1:
        # primary reactions only
        rcella = -(K2+K3)*cella    # active cellulose rate
        rgas = K3*cella            # gas rate
        rtar = K2*cella            # tar rate
        rchar = K3*cella           # char rate
        ncella = cella + rcella*dt      # update active cellulose concentration
        ngas = gas + rgas*(1-x)*dt      # update gas concentration
        ntar = tar + rtar*dt            # update tar concentration
        nchar = char + rchar*x*dt       # update char concentration
    elif s == 2:
        # primary and secondary reactions
        rcella = -(K2+K3)*cella    # active cellulose rate
        rgas = K3*cella + K4*tar   # gas rate
        rtar = K2*cella - K4*tar   # tar rate
        rchar = K3*cella           # char rate
        ncella = cella + rcella*dt      # update active cellulose concentration
        ngas = gas + rgas*(1-x)*dt      # update gas concentration
        ntar = tar + rtar*dt            # update tar concentration
        nchar = char + rchar*x*dt       # update char concentration

    # return new cellulose, active cellulose, gas, tar, char mass concentrations
    return ncella, ngas, ntar, nchar


def millerhemi_noR1(hemia, gas, tar, char, T, dt, s=1):
    """
    Hemicellulose reactions from Miller and Bellan 1997 paper for biomass
    pyrolysis.

    Parameters
    ----------
    hemia = active hemicellulose concentration, kg/m^3
    gas = gas concentration, kg/m^3
    tar = tar concentration, kg/m^3
    char = char concentration, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nhemia = new active hemicellulose concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A2 = 8.75e15;   E2 = 202.4    # hemia -> tar
    A3 = 2.6e11;    E3 = 145.7    # hemia -> x*char + (1-x)*gas
    A4 = 4.28e6;    E4 = 108      # tar -> gas
    R = 0.008314    # universal gas constant, kJ/mol*K
    x = 0.60        # char formation mass ratio for hemicellulose

    # reaction rate constant for each reaction, 1/s
    K2 = A2 * np.exp(-E2 / (R * T))  # hemia -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # hemia -> x*char + (1-x)*gas
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas

    if s == 1:
        # primary reactions only
        rhemia = -(K2+K3)*hemia    # active hemicellulose rate
        rgas = K3*hemia            # gas rate
        rtar = K2*hemia            # tar rate
        rchar = K3*hemia           # char rate
        nhemia = hemia + rhemia*dt      # update active hemicellulose concentration
        ngas = gas + rgas*(1-x)*dt      # update gas concentration
        ntar = tar + rtar*dt            # update tar concentration
        nchar = char + rchar*x*dt       # update char concentration
    elif s == 2:
        # primary and secondary reactions
        rhemia = -(K2+K3)*hemia    # active hemicellulose rate
        rgas = K3*hemia + K4*tar   # gas rate
        rtar = K2*hemia - K4*tar   # tar rate
        rchar = K3*hemia           # char rate
        nhemia = hemia + rhemia*dt      # update active hemicellulose concentration
        ngas = gas + rgas*(1-x)*dt      # update gas concentration
        ntar = tar + rtar*dt            # update tar concentration
        nchar = char + rchar*x*dt       # update char concentration

    # return new hemicellulose, active hemicellulose, gas, tar, char mass concentrations
    return nhemia, ngas, ntar, nchar


def millerlig_noR1(liga, gas, tar, char, T, dt, s=1):
    """
    Lignin reactions from Miller and Bellan 1997 paper for biomass pyrolysis.

    Parameters
    ----------
    liga = active lignin concentration, kg/m^3
    gas = gas concentration, kg/m^3
    tar = tar concentration, kg/m^3
    char = char concentration, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nliga = new active lignin concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A2 = 8.75e15;   E2 = 202.4    # liga -> tar
    A3 = 2.6e11;    E3 = 145.7    # liga -> x*char + (1-x)*gas
    A4 = 4.28e6;    E4 = 108      # tar -> gas
    R = 0.008314    # universal gas constant, kJ/mol*K
    x = 0.75        # char formation mass ratio for lignin

    # reaction rate constant for each reaction, 1/s
    K2 = A2 * np.exp(-E2 / (R * T))  # liga -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # liga -> x*char + (1-x)*gas
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas

    if s == 1:
        # primary reactions only
        rliga = -(K2+K3)*liga   # active lignin rate
        rgas = K3*liga          # gas rate
        rtar = K2*liga          # tar rate
        rchar = K3*liga         # char rate
        nliga = liga + rliga*dt     # update active lignin concentration
        ngas = gas + rgas*(1-x)*dt  # update gas concentration
        ntar = tar + rtar*dt        # update tar concentration
        nchar = char + rchar*x*dt   # update char concentration
    elif s == 2:
        # primary and secondary reactions
        rliga = -(K2+K3)*liga   # active lignin rate
        rgas = K3*liga + K4*tar # gas rate
        rtar = K2*liga - K4*tar # tar rate
        rchar = K3*liga         # char rate
        nliga = liga + rliga*dt     # update active lignin concentration
        ngas = gas + rgas*(1-x)*dt  # update gas concentration
        ntar = tar + rtar*dt        # update tar concentration
        nchar = char + rchar*x*dt   # update char concentration

    # return vector for each concentration, kg/m^3
    return nliga, ngas, ntar, nchar

# Products from Cellulose, Hemicellulose, and Lignin
# ------------------------------------------------------------------------------

cella = np.ones((2, nt))*wtcell
hemia = np.ones((2, nt))*wthemi
liga = np.ones((2, nt))*wtlig
gas = np.zeros((6, nt))
tar = np.zeros((6, nt))
char = np.zeros((6, nt))

for i in range(1, nt):
    cella[0, i], gas[0, i], tar[0, i], char[0, i] = millercell_noR1(cella[0, i-1], gas[0, i-1], tar[0, i-1], char[0, i-1], T, dt)
    cella[1, i], gas[1, i], tar[1, i], char[1, i] = millercell_noR1(cella[1, i-1], gas[1, i-1], tar[1, i-1], char[1, i-1], T, dt, s=2)
    hemia[0, i], gas[2, i], tar[2, i], char[2, i] = millerhemi_noR1(hemia[0, i-1], gas[2, i-1], tar[2, i-1], char[2, i-1], T, dt)
    hemia[1, i], gas[3, i], tar[3, i], char[3, i] = millerhemi_noR1(hemia[1, i-1], gas[3, i-1], tar[3, i-1], char[3, i-1], T, dt, s=2)
    liga[0, i], gas[4, i], tar[4, i], char[4, i] = millerlig_noR1(liga[0, i-1], gas[4, i-1], tar[4, i-1], char[4, i-1], T, dt)
    liga[1, i], gas[5, i], tar[5, i], char[5, i] = millerlig_noR1(liga[1, i-1], gas[5, i-1], tar[5, i-1], char[5, i-1], T, dt, s=2)

# sum products into groups of wood, gas, tar, and char
wood1 = cella[0] + hemia[0] + liga[0]
wood2 = cella[1] + hemia[1] + liga[1]

gas1 = gas[0] + gas[2] + gas[4]
gas2 = gas[1] + gas[3] + gas[5]

tar1 = tar[0] + tar[2] + tar[4]
tar2 = tar[1] + tar[3] + tar[5]

char1 = char[0] + char[2] + char[4]
char2 = char[1] + char[3] + char[5]

# Plot Results
# ------------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, wood1, lw=2, label='wood1')
py.plot(t, gas1, lw=2, label='gas1')
py.plot(t, tar1, lw=2, label='tar1')
py.plot(t, char1, lw=2, label='char1')
py.title('Miller 1997 primary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Conversion (mass fraction)')
py.legend(loc='best', numpoints=1, fontsize=11, frameon=False)
py.grid()

py.figure(2)
py.plot(t, wood2, lw=2, label='wood2')
py.plot(t, gas2, lw=2, label='gas2')
py.plot(t, tar2, lw=2, label='tar2')
py.plot(t, char2, lw=2, label='char2')
py.title('Miller 1997 primary and secondary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Conversion (mass fraction)')
py.legend(loc='best', numpoints=1, fontsize=11, frameon=False)
py.grid()
