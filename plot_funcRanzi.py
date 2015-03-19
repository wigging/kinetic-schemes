"""
Plot the kinetic reactions of biomass pyrolysis for the Ranzi kinetic scheme.
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
import funcRanzi as kn

py.close('all')

# Parameters from Papadikis 2010a
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.01   # time step, delta t
tmax = 4   # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
p = len(t)  # total number of time steps

# Calculate Concentrations of Chemical Species
#------------------------------------------------------------------------------

# vectors for wood, gas, tar, char concentrations as a density, kg/m^3
pw = np.zeros(len(t))   # wood 
pg = np.zeros(len(t))   # gas
pt = np.zeros(len(t))   # tar
pc = np.zeros(len(t))   # char

pw[:] = rhow    # initial wood concentration as density

# array of chemical species concentrations as a density, kg/m^3
pcell, pcella, pg1, pg2, plvg = kn.cell2(Tinf, pw, dt, p)

# concentration as percent relative to original wood, %
cell = pcell/rhow*100       # CELL
cella = pcella/rhow*100     # CELLA
g1 = pg1/rhow*100           # G1
g2 = pg2/rhow*100           # G2
lvg = plvg/rhow*100         # LVG

total = cell + cella + g1 + g2 + lvg # mass balance for CELL, CELLA, G1, G2, LVG

h2o = 5*g1/(5+6)            # H2O
char = 6*g1/(5+6)           # Char

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(1)
py.plot(t, cell, label='cell')
py.plot(t, cella, label='cella')
py.plot(t, g1, label='g1')
py.plot(t, g2, label='g2')
py.plot(t, lvg, label='lvg')
py.plot(t, total, label='total')
py.plot(t, h2o, label='h2o')
py.plot(t, char, label='char')
py.title('Cellulose Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Conversion (% dry basis)')
py.legend(loc='best', numpoints=1)
py.show()
