"""
Simple example of mass kinetic scheme to mole kinetic scheme.
"""

import numpy as np
import matplotlib.pyplot as py

# Parameters
# ------------------------------------------------------------------------------

rho = 540   # density of pine, kg/m^3
T = 773     # temperature for rate constants, K

dt = 0.01                               # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Function
# ------------------------------------------------------------------------------

def cpc(wood, tar, gas, T, dt):
    """
    Calculate conversion of wood to tar and gas using Chan/Liden kinetics.
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A2 = 2.0e8;     E2 = 133    # wood -> tar   from Chan 1985
    A4 = 4.28e6;    E4 = 107.5  # tar -> gas    from Liden 1988
    R = 0.008314    # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K2 = A2 * np.exp(-E2 / (R * T))   # wood -> tar
    K4 = A4 * np.exp(-E4 / (R * T))   # tar -> gas

    # primary and secondary reactions
    rw = -K2*wood           # wood rate
    rt = K2*wood - K4*tar   # tar rate
    rg = K4*tar             # gas rate
    nwood = wood + rw*dt            # update wood concentration
    ntar = tar + rt*dt              # update tar concentration
    ngas = gas + rg*dt              # update gas concentration

    # return new mass concentrations of products, kg/m^3
    return nwood, ntar, ngas

# Mass based approach with C = kg/m^3 
# ------------------------------------------------------------------------------

# store concentrations from primary and secondary reactions at each time step
# concentrations calculated on a mass basis such as kg/m^3
wood = np.ones(nt)*rho    # wood concentration vector
tar = np.zeros(nt)        # tar concentration vector
gas = np.zeros(nt)        # gas concentration vector

# calculate products at each time step
for i in range(1, nt):
    wood[i], tar[i], gas[i] = cpc(wood[i-1], tar[i-1], gas[i-1], T, dt)

# mass balance to check total mass concentration at each time step
mass = wood + tar + gas
print('total mass is', mass)

# Mole based approach with C = mol/m^3
# ------------------------------------------------------------------------------

# assume molecular weight of wood is 50 g/mol
mw = 50/1000    # kg/mol

# store concentrations from primary and secondary reactions at each time step
# concentrations calculated on a mass basis such as kg/m^3
wood2 = np.ones(nt)*rho/mw    # wood concentration vector
tar2 = np.zeros(nt)           # tar concentration vector
gas2 = np.zeros(nt)           # gas concentration vector

# calculate products at each time step
for i in range(1, nt):
    wood2[i], tar2[i], gas2[i] = cpc(wood2[i-1], tar2[i-1], gas2[i-1], T, dt)

# mole balance to check total mole concentration at each time step
moles = wood2 + tar2 + gas2
print('total moles is', moles)

# Plot
# ------------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, wood, lw=2, label='wood')
py.plot(t, tar, lw=2, label='tar')
py.plot(t, gas, lw=2, label='gas')
py.title('Mass basis at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration (kg/m^3)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t, wood2, lw=2, label='wood')
py.plot(t, tar2, lw=2, label='tar')
py.plot(t, gas2, lw=2, label='gas')
py.title('Mole basis at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration (mol/m^3)')
py.legend(loc='best', numpoints=1)
py.grid()

