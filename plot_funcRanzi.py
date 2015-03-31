"""
Plot the kinetic reactions of biomass pyrolysis for the Ranzi kinetic scheme.
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
import funcRanzi as kn

py.close('all')

# Parameters
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# weight percent (%) cellulose, hemicellulose, lignin for beech wood
wtcell  = 48
wthemi = 28
wtlig = 24

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.01                               # time step, delta t
tmax = 4                                # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Calculate Products from Cellulose Reactions
#------------------------------------------------------------------------------

# array of product concentrations as a density, kg/m^3
pcell = kn.cell(rhow, wtcell, Tinf, dt, nt)

# convert product concentrations to mass fraction, (-)    
cell = pcell[0]/rhow    # CELL
cella = pcell[1]/rhow   # CELLA
g1 = pcell[2]/rhow      # G1
g2 = pcell[3]/rhow      # G2
lvg = pcell[4]/rhow     # LVG

total = cell + cella + g1 + g2 + lvg # mass balance for CELL, CELLA, G1, G2, LVG

h2o = 5*g1/(5+6)            # H2O
char = 6*g1/(5+6)           # Char

total2 = cell+cella+h2o+char+g2+lvg # check mass balance again

# Calculate Products from Hemicellulose Reactions
#------------------------------------------------------------------------------

# array of product concentrations as a density, kg/m^3
phce = kn.hce(rhow, wthemi, Tinf, dt, nt)

# convert product concentrations to mass fraction, (-)    
hce = phce[0]/rhow
g11 = phce[1]/rhow
g22 = phce[2]/rhow
g33 = phce[3]/rhow
g44 = phce[4]/rhow
xyl = phce[5]/rhow

# Calculate Products from Lignin Reactions
#------------------------------------------------------------------------------

# array of product concentrations as a density, kg/m^3
pligc = kn.ligc(rhow, wtlig, Tinf, dt, nt)

# convert product concentrations to mass fraction, (-)
ligc = pligc[0]/rhow
g111 = pligc[1]/rhow
g222 = pligc[2]/rhow

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
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(2)
py.plot(t, hce, label='hce')
py.plot(t, g11, label='g1')
py.plot(t, g22, label='g2')
py.plot(t, g33, label='g3')
py.plot(t, g44, label='g4')
py.plot(t, xyl, label='xyl')
py.title('Hemicellulose Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(3)
py.plot(t, ligc, label='ligc')
py.plot(t, g111, label='g1')
py.plot(t, g222, label='g2')
py.title('Lignin Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.show()
