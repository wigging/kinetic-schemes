"""
Plot gas, tar, char from primary and secondary reactions as determined by the
Janse 2000 kinetic scheme for biomass pyrolysis.

Reference:
Janse, Westerhout, Prins, 2000. Chem. Eng. Process., 39, pp 239-252.
"""

import numpy as np
import matplotlib.pyplot as py

# Parameters
# ------------------------------------------------------------------------------

T = 773  # temperature for rate constants, K

dt = 0.01                               # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Function for Janse 2000 Kinetic Scheme
# ------------------------------------------------------------------------------

def janse(wood, gas, tar, char, T, dt, s=1):
    """
    Primary and secondary kinetic reactions from Table 1 in Janse 2000 paper.
    Also see paper by Liden 1988.

    Parameters
    ----------
    wood = wood concentration, kg/m^3
    gas = gas concentation, kg/m^3
    tar = tar concentation, kg/m^3
    char = char concentation, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nwood = new wood concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.11e11;   E1 = 177     # wood -> gas
    A2 = 9.28e9;    E2 = 149     # wood -> tar
    A3 = 3.05e7;    E3 = 125     # wood -> char
    A4 = 8.6e4;     E4 = 87.8    # tar -> gas
    A5 = 7.7e4;     E5 = 87.8    # tar -> char
    R = 0.008314    # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas
    K5 = A5 * np.exp(-E5 / (R * T))  # tar -> char

    if s == 1:
        # primary reactions only
        rw = -(K1+K2+K3)*wood     # wood rate
        rg = K1*wood              # gas rate
        rt = K2*wood              # tar rate
        rc = K3*wood              # char rate
        nwood = wood + rw*dt     # update wood concentration
        ngas = gas + rg*dt     # update gas concentration
        ntar = tar + rt*dt     # update tar concentration
        nchar = char + rc*dt     # update char concentration
    elif s == 2:
        # primary and secondary reactions
        rw = -(K1+K2+K3)*wood         # wood rate
        rg = K1*wood + K4*tar          # gas rate
        rt = K2*wood - (K4+K5)*tar     # tar rate
        rc = K3*wood + K5*tar          # char rate
        nwood = wood + rw*dt        # update wood concentration
        ngas = gas + rg*dt        # update gas concentration
        ntar = tar + rt*dt        # update tar concentration
        nchar = char + rc*dt        # update char concentration

    # return new wood, gas, tar, char mass concentrations, kg/m^3
    return nwood, ngas, ntar, nchar

# Products from Kinetic Scheme
# ------------------------------------------------------------------------------

# store concentrations from primary reactions at each time step
# concentrations reported on a mass basis as kg/m^2
wood = np.ones(nt)  # wood concentration vector
gas = np.zeros(nt)  # gas concentration vector
tar = np.zeros(nt)  # tar concentration vector
char = np.zeros(nt) # char concentration vector

# store concentrations from primary and secondary reactions at each time step
wood2 = np.ones(nt)     # wood concentration vector
gas2 = np.zeros(nt)     # gas concentration vector
tar2 = np.zeros(nt)     # tar concentration vector
char2 = np.zeros(nt)    # char concentration vector

# products from primary reactions only
for i in range(1, nt):
    wood[i], gas[i], tar[i], char[i] = janse(wood[i-1], gas[i-1], tar[i-1], char[i-1], T, dt)

# products from primary and secondary reactions
for i in range(1, nt):
    wood2[i], gas2[i], tar2[i], char2[i] = janse(wood2[i-1], gas2[i-1], tar2[i-1], char2[i-1], T, dt, s=2)

# Plot Results
# ------------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, wood, lw=2, label='wood')
py.plot(t, gas, lw=2, label='gas')
py.plot(t, tar, lw=2, label='tar')
py.plot(t, char, lw=2, label='char')
py.title('Janse 2000 primary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t, wood2, lw=2, label='wood')
py.plot(t, gas2, lw=2, label='gas')
py.plot(t, tar2, lw=2, label='tar')
py.plot(t, char2, lw=2, label='char')
py.title('Janse 2000 primary and secondary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.legend(loc='best', numpoints=1)
py.grid()
