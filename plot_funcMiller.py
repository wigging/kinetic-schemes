"""
Plot the kinetic reactions of biomass pyrolysis for the Miller kinetic scheme.
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
import funcMiller as kn

py.close('all')

# Parameters from Papadikis 2010a
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.0001   # time step, delta t
tmax = 25   # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
p = len(t)  # total number of time steps

# Calculate Concentrations of Chemical Species
#------------------------------------------------------------------------------

# vector for initial wood concentrations as a density, kg/m^3
pw = np.ones(p)*rhow

# array of product concentrations as a density, kg/m^3
pcell, pcella, pgas1, ptar1, pchar1 = kn.cell(Tinf, pw, dt, p)
phemi, phemia, pgas2, ptar2, pchar2 = kn.hemi(Tinf, pw, dt, p)
plig, pliga, pgas3, ptar3, pchar3 = kn.lig(Tinf, pw, dt, p)

# sum components to wood, gas, tar, char groups, kg/m^3
pwood = pcell + phemi + plig
pactive = pcella + phemia + pliga
pgas = pgas1 + pgas2 + pgas3
ptar = ptar1 + ptar2 + ptar3
pchar = pchar1 + pchar2 + pchar3

# concentration as percent relative to original wood, %
wood = pwood/rhow*100       # wood
active = pactive/rhow*100   # active
gas = pgas/rhow*100         # gas
tar = ptar/rhow*100         # tar
char = pchar/rhow*100       # char

total = wood + active + gas + tar + char # mass balance

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(1)
py.plot(t, wood, label='wood')
py.plot(t, gas, label='gas')
py.plot(t, tar, label='tar')
py.plot(t, char, label='char')
#py.plot(t, total, label='total')
py.title('Cell, Hemi, Lig Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Conversion (% dry basis)')
py.legend(loc='best', numpoints=1)
py.show()
