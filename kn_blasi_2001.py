"""
Plot gas, tar, char from primary reactions as determined by the Blasi Branca 2001
kinetic scheme for biomass pyrolysis.

Reference:
Blasi, Branca, 2001. Ind. Eng. Chem. Res., 40, pp 5547-5556.
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

# Function for Blasi Branca 2001 Kinetic Scheme
# ------------------------------------------------------------------------------

def blasibranca(pw, pg, pt, pc, T, dt):
    """
    Primary kinetic reactions from Table 1 in Blasi and Branca 2001 paper.

    Parameters
    ----------
    pw = wood concentration, kg/m^3
    pg = gas concentation, kg/m^3
    pt = tar concentation, kg/m^3
    pc = char concentation, kg/m^3
    T = temperature, K
    dt = time step, s

    Returns
    -------
    nw = new wood concentration, kg/m^3
    ng = new gas concentration, kg/m^3
    nt = new tar concentration, kg/m^3
    nc = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 4.38e9;    E1 = 152.7    # wood -> gas
    A2 = 1.08e10;   E2 = 148      # wood -> tar
    A3 = 3.27e6;    E3 = 111.7    # wood -> char
    R = 0.008314 # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char

    # primary reactions only
    rw = -(K1+K2+K3)*pw # wood rate
    rg = K1*pw          # gas rate
    rt = K2*pw          # tar rate
    rc = K3*pw          # char rate
    nw = pw + rw*dt # update wood concentration
    ng = pg + rg*dt # update gas concentration
    nt = pt + rt*dt # update tar concentration
    nc = pc + rc*dt # update char concentration

    # return new wood, gas, tar, char as mass concentrations, kg/m^3
    return nw, ng, nt, nc

# Products from Kinetic Scheme
# ------------------------------------------------------------------------------

# store concentrations from primary reactions at each time step
# concentrations reported on a mass basis as kg/m^2
wood = np.ones(nt)  # wood concentration vector
gas = np.zeros(nt)  # gas concentration vector
tar = np.zeros(nt)  # tar concentration vector
char = np.zeros(nt) # char concentration vector

# products from primary reactions only
for i in range(1, nt):
    wood[i], gas[i], tar[i], char[i] = blasibranca(wood[i-1], gas[i-1], tar[i-1], char[i-1], T, dt)

# Plot Results
# ------------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, wood, lw=2, label='wood')
py.plot(t, gas, lw=2, label='gas')
py.plot(t, tar, lw=2, label='tar')
py.plot(t, char, lw=2, label='char')
py.title('Blasi 2001 primary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration (mass fraction)')
py.legend(loc='best', numpoints=1)
py.grid()
