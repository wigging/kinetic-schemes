"""
Plot gas, tar, char from primary reactions as determined by the Thurner 1981
kinetic scheme for biomass pyrolysis.

Reference:
Thurner, Mann, 1981. Ind. Eng. Chem. Process Des. Dev., 20, pp 482-488.
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

# Function for Thurner 1981 Kinetic scheme
# ------------------------------------------------------------------------------

def thurner(wood, gas, tar, char, T, dt):
    """
    Primary kinetic reactions from Thurner 1981 paper for biomass pyrolysis.

    Parameters
    ----------
    wood = wood concentration, kg/m^3
    gas = gas concentation, kg/m^3
    tar = tar concentation, kg/m^3
    char = char concentation, kg/m^3
    T = temperature, K
    dt = time step, s

    Returns
    -------
    nwood = new wood concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.44e4;    E1 = 88.6     # wood -> gas
    A2 = 4.13e6;    E2 = 112.7    # wood -> tar
    A3 = 7.38e5;    E3 = 106.5    # wood -> char
    R = 0.008314    # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char

    # primary reactions only
    rw = -(K1+K2+K3) * wood     # wood rate
    rg = K1 * wood              # gas rate
    rt = K2 * wood              # tar rate
    rc = K3 * wood              # char rate
    nwood = wood + rw*dt    # update wood concentration
    ngas = gas + rg*dt      # update gas concentration
    ntar = tar + rt*dt      # update tar concentration
    nchar = char + rc*dt    # update char concentration

    # return new wood, char, gas, tar as mass concentrations, kg/m^3
    return nwood, ngas, ntar, nchar

# Products from Kinetic Scheme
# ------------------------------------------------------------------------------

# concentrations reported on a mass basis in kg/m^3
# pw = wood concentration, pg = gas concentration, pt = tar concentration,
# pc = char concentration

# store concentrations from primary reactions at each time step
wood = np.ones(nt)      # wood concentration vector
gas = np.zeros(nt)      # gas concentration vector
tar = np.zeros(nt)      # tar concentration vector
char = np.zeros(nt)     # char concentration vector

# products from primary reactions only
for i in range(1, nt):
    wood[i], gas[i], tar[i], char[i] = thurner(wood[i-1], gas[i-1], tar[i-1], char[i-1], T, dt)

# Plot Results
# ------------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, wood, lw=2, label='wood')
py.plot(t, gas, lw=2, label='gas')
py.plot(t, tar, lw=2, label='tar')
py.plot(t, char, lw=2, label='char')
py.title('Thurner 1981 primary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration (mass fraction)')
py.legend(loc='best', numpoints=1)
py.grid()
