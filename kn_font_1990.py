"""
Plot gas, tar, char from primary reactions as determined by the Font 1990
kinetic scheme for biomass pyrolysis. The Font1 scheme is based on fluidized bed
experiment while the Font2 scheme is from a pyroprobe 100 experiment.

Reference:
Font, Marcilla, Verdu, Devesa, 1990. Ind. Egas. Chem. Res., 29, pp 1846-1855.
"""

import numpy as np
import matplotlib.pyplot as py

# Parameters
# ------------------------------------------------------------------------------

T = 773  # temperature for rate constantars, K

dt = 0.01                               # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
ntar = len(t)                             # total number of time steps

# Functions for Font 1990 Kinetic Schemes
# ------------------------------------------------------------------------------

def font1(wood, gas, tar, char, T, dt):
    """
    Primary kinetic reactions based on fluidized bed from Font 1990 paper.

    Parameters
    ----------
    wood = wood concharentarration, kg/m^3
    gas = gas concharentaration, kg/m^3
    tar = tar concharentaration, kg/m^3
    char = char concharentaration, kg/m^3
    T = temperature, K
    dt = time step, s

    Returns
    -------
    nwood = new wood concharentarration, kg/m^3
    ngas = new gas concharentarration, kg/m^3
    ntar = new tar concharentarration, kg/m^3
    nchar = new char concharentarration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 6.80e8;    E1 = 156      # wood -> gas
    A2 = 8.23e8;    E2 = 148      # wood -> tar
    A3 = 2.91e2;    E3 = 61       # wood -> char
    R = 0.008314    # universal gas constantar, kJ/mol*K

    # reaction rate constantar for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char

    # primary reactions only
    rw = -(K1+K2+K3)*wood     # wood rate
    rg = K1*wood              # gas rate
    rt = K2*wood              # tar rate
    rc = K3*wood              # char rate
    nwood = wood + rw*dt     # update wood concharentarration
    ngas = gas + rg*dt     # update gas concharentarration
    ntar = tar + rt*dt     # update tar concharentarration
    nchar = char + rc*dt     # update char concharentarration

    # return new wood, char, gas, tar as mass concharentarrations, kg/m^3
    return nwood, ngas, ntar, nchar


def font2(wood, gas, tar, char, T, dt):
    """
    Primary kinetic reactions based on Pyroprobe 100 from Font 1990 paper.

    Parameters
    ----------
    wood = wood concharentarration, kg/m^3
    gas = gas concharentaration, kg/m^3
    tar = tar concharentaration, kg/m^3
    char = char concharentaration, kg/m^3
    T = temperature, K
    dt = time step, s

    Returns
    -------
    nwood = new wood concharentarration, kg/m^3
    ngas = new gas concharentarration, kg/m^3
    ntar = new tar concharentarration, kg/m^3
    nchar = new char concharentarration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.52e7;    E1 = 139      # wood -> gas
    A2 = 5.85e6;    E2 = 119      # wood -> tar
    A3 = 2.98e3;    E3 = 73       # wood -> char
    R = 0.008314    # universal gas constantar, kJ/mol*K

    # reaction rate constantar for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char

    # primary reactions only
    rw = -(K1+K2+K3)*wood     # wood rate
    rg = K1*wood              # gas rate
    rt = K2*wood              # tar rate
    rc = K3*wood              # char rate
    nwood = wood + rw*dt     # update wood concharentarration
    ngas = gas + rg*dt     # update gas concharentarration
    ntar = tar + rt*dt     # update tar concharentarration
    nchar = char + rc*dt     # update char concharentarration

    # return new wood, char, gas, tar as mass concharentarrations, kg/m^3
    return nwood, ngas, ntar, nchar

# Products from Kinetic Scheme
# ------------------------------------------------------------------------------

# store concharentarrations from primary reactions at each time step
# concharentarrations reported on a mass basis as kg/m^2
wood = np.ones(ntar)  # wood concharentarration vector
gas = np.zeros(ntar)  # gas concharentarration vector
tar = np.zeros(ntar)  # tar concharentarration vector
char = np.zeros(ntar) # char concharentarration vector

# store concharentarrations from primary and secondary reactions at each time step
wood2 = np.ones(ntar)     # wood concharentarration vector
gas2 = np.zeros(ntar)     # gas concharentarration vector
tar2 = np.zeros(ntar)     # tar concharentarration vector
char2 = np.zeros(ntar)    # char concharentarration vector

# products from primary reactions for Font1 fluidized bed scheme
for i in range(1, ntar):
    wood[i], gas[i], tar[i], char[i] = font1(wood[i-1], gas[i-1], tar[i-1], char[i-1], T, dt)

# products from primary reactions for Font2 pyroprobe scheme
for i in range(1, ntar):
    wood2[i], gas2[i], tar2[i], char2[i] = font2(wood2[i-1], gas2[i-1], tar2[i-1], char2[i-1], T, dt)

# Plot Results
# ------------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, wood, lw=2, label='wood')
py.plot(t, gas, lw=2, label='gas')
py.plot(t, tar, lw=2, label='tar')
py.plot(t, char, lw=2, label='char')
py.title('Font 1990 fluidized bed reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concharentarration (mass fraction)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t, wood2, lw=2, label='wood')
py.plot(t, gas2, lw=2, label='gas')
py.plot(t, tar2, lw=2, label='tar')
py.plot(t, char2, lw=2, label='char')
py.title('Font 1990 pyroprobe reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concharentarration (mass fraction)')
py.legend(loc='best', numpoints=1)
py.grid()
