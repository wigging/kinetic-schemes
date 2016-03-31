"""
Plot gas, tar, char from primary and secondary reactions as determined by the
Chan 1985 kinetic scheme for biomass pyrolysis. Note that secondary reactions do
not reduce tar, it produces more tar labeled as tar2.

Reference:
Chan, Kelbon, Krieger, 1985. Fuel, 64, pp 1505–1513.
"""

import numpy as np
import matplotlib.pyplot as py

# Parameters
# ------------------------------------------------------------------------------

T = 773     # temperature for rate constants, K
mc = 0.10   # moisture content as fraction, (-)

dt = 0.01                               # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Function for Chan 1985 Kinetic Scheme
# ------------------------------------------------------------------------------

def chan(wood, gas, tar, char, water, vapor, T, dt, s=1):
    """
    Primary and secondary reactions with moisture reaction from Table 2 in the
    Chan 1985 paper.

    Parameters
    ----------
    wood = wood concentration, kg/m^3
    gas = gas concentation, kg/m^3
    tar = tar concentation, kg/m^3
    char = char concentation, kg/m^3
    water = moisture content concentration, kg/m^3
    vapor = water vapor concentration, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nwood = new wood concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    nwater = new moisture concentration, kg/m^3
    nvapor = new water vapor concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.3e8;     E1 = 140    # wood -> gas1
    A2 = 2.0e8;     E2 = 133    # wood -> tar1
    A3 = 1.08e7;    E3 = 121    # wood -> char
    A4 = 5.13e6;    E4 = 87.9   # moisture -> water vapor
    A5 = 1.48e6;    E5 = 144    # tar1 -> gas2 + tar2
    R = 0.008314    # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas1
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar1
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T))  # moisture -> water vapor
    K5 = A5 * np.exp(-E5 / (R * T))  # tar1 -> gas2 + tar2

    if s == 1:
        # primary reactions only
        rw = -(K1+K2+K3)*wood     # wood rate
        rg = K1*wood              # gas rate
        rt = K2*wood              # tar rate
        rc = K3*wood              # char rate
        rm = -K4*water             # moisture content rate
        rv = K4*water              # water vapor rate
        nwood = wood + rw*dt     # update wood concentration
        ngas = gas + rg*dt     # update gas concentration
        ntar = tar + rt*dt     # update tar concentration
        nchar = char + rc*dt     # update char concentration
        nwater = water + rm*dt     # update moisture content concentration
        nvapor = vapor + rv*dt     # update water vapor concentration
    elif s == 2:
        # primary and secondary reactions
        rw = -(K1+K2+K3)*wood     # wood rate
        rg = K1*wood + K5*tar      # gas rate
        rt = K2*wood - K5*tar      # tar rate
        rc = K3*wood              # char rate
        rm = -K4*water             # moisture content rate
        rv = K4*water              # water vapor rate
        nwood = wood + rw*dt     # update wood concentration
        ngas = gas + rg*dt     # update gas concentration
        ntar = tar + rt*dt     # update tar concentration
        nchar = char + rc*dt     # update char concentration
        nwater = water + rm*dt     # update moisture content concentration
        nvapor = vapor + rv*dt     # update water vapor concentration

    # return wood, gas, tar, char, moisture, water vapor as mass concentrations
    return nwood, ngas, ntar, nchar, nwater, nvapor

# Products from Kinetic Scheme
# ------------------------------------------------------------------------------

# store concentrations from primary reactions at each time step
# concentrations reported on a mass basis as kg/m^2
wood = np.ones(nt)      # wood concentration vector
gas = np.zeros(nt)      # gas concentration vector
tar = np.zeros(nt)      # tar concentration vector
char = np.zeros(nt)     # char concentration vector
water = np.ones(nt)*mc     # moisture content concentration vector
vapor = np.zeros(nt)       # water vapor concentration vector

# store concentrations from primary and secondary reactions at each time step
wood2 = np.ones(nt)     # wood concentration vector
gas2 = np.zeros(nt)     # gas concentration vector
tar2 = np.zeros(nt)     # tar concentration vector
char2 = np.zeros(nt)    # char concentration vector
water2 = np.ones(nt)*mc    # moisture content concentration vector
vapor2 = np.zeros(nt)      # water vapor concentration vector

# products from primary reactions only
for i in range(1, nt):
    wood[i], gas[i], tar[i], char[i], water[i], vapor[i] = chan(wood[i-1], gas[i-1], tar[i-1], char[i-1], water[i-1], vapor[i-1], T, dt)

# products from primary and secondary reactions
for i in range(1, nt):
    wood2[i], gas2[i], tar2[i], char2[i], water2[i], vapor2[i] = chan(wood2[i-1], gas2[i-1], tar2[i-1], char2[i-1], water2[i-1], vapor2[i-1], T, dt, s=2)

# Plot Results
# ------------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, wood, lw=2, label='wood')
py.plot(t, gas, lw=2, label='gas')
py.plot(t, tar, lw=2, label='tar')
py.plot(t, char, lw=2, label='char')
py.plot(t, water, lw=2, label='moisture')
py.plot(t, vapor, lw=2, label='water vapor')
py.title('Chan 1985 primary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration (mass fraction)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t, wood2, lw=2, label='wood')
py.plot(t, gas2, lw=2, label='gas')
py.plot(t, tar2, lw=2, label='tar')
py.plot(t, char2, lw=2, label='char')
py.plot(t, water2, lw=2, label='moisture')
py.plot(t, vapor2, lw=2, label='water vapor')
py.title('Chan 1985 primary and secondary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration (mass fraction)')
py.legend(loc='best', numpoints=1)
py.grid()
