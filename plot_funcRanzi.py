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

# Products from Cellulose, Hemicellulose, Lignin Reactions
#------------------------------------------------------------------------------

# array of product concentrations as a density, kg/m^3
pcell = kn.cell(rhow, wtcell, Tinf, dt, nt)

# products as mass fraction, (-)    
cell = pcell[0]/rhow    # CELL
cella = pcell[1]/rhow   # CELLA
g1 = pcell[2]/rhow      # G1
g4 = pcell[3]/rhow      # G4
lvg = pcell[4]/rhow     # LVG
h2o = pcell[5]/rhow     # H2O
char = pcell[6]/rhow    # Char
haa = pcell[7]/rhow     # HAA
glyox = pcell[8]/rhow   # Glyox
c2h4o = pcell[9]/rhow   # C2H4O
hmfu = pcell[10]/rhow   # HMFU
c3h6o = pcell[11]/rhow  # C3H6O
co2 = pcell[12]/rhow    # CO2
h2 = pcell[13]/rhow     # H2
ch2o = pcell[14]/rhow   # CH2O
co = pcell[15]/rhow     # CO
ch4 = pcell[16]/rhow    # CH4
hcooh = pcell[17]/rhow  # HCOOH

# mass balances
total = cell + cella + g1 + g4 + lvg # for CELL, CELLA, G1, G$, LVG

total2 = cell + cella + lvg + h2o + char + haa + glyox + c2h4o + hmfu + c3h6o + co2 + \
         h2 + ch2o + co + ch4 + hcooh  # approach 2

# Products from Hemicellulose Reactions
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

# Products from Lignin Reactions
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
py.plot(t, g4, label='g4')
py.plot(t, lvg, label='lvg')
py.plot(t, total, label='total')
py.title('Cellulose Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(2)
py.plot(t, h2o, label='$\mathrm{H_{2}O}$')
py.plot(t, char, label='char')
py.plot(t, haa, label='haa')
py.plot(t, glyox, label='glyox')
py.plot(t, c2h4o, label='$\mathrm{C_{2}H_{4}O}$')
py.plot(t, hmfu, label='HMFU')
py.plot(t, c3h6o, label='$\mathrm{C_{3}H_{6}O}$')
py.plot(t, co2, ls='--', label='$\mathrm{CO_2}$')
py.plot(t, h2, ls='--', label='$\mathrm{H_2}$')
py.plot(t, ch2o, ls='--', label='$\mathrm{CH_2O}$')
py.plot(t, co, ls='--', label='$\mathrm{CO}$')
py.plot(t, ch4, ls='--', label='$\mathrm{CH_4}$')
py.plot(t, hcooh, ls='--', label='$\mathrm{HCOOH}$')
py.plot(t, lvg, ls='--', label='lvg')
py.plot(t, total2, ls='--', label='total2')
py.title('Cellulose Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

#py.figure(2)
#py.plot(t, hce, label='hce')
#py.plot(t, g11, label='g1')
#py.plot(t, g22, label='g2')
#py.plot(t, g33, label='g3')
#py.plot(t, g44, label='g4')
#py.plot(t, xyl, label='xyl')
#py.title('Hemicellulose Reactions at T = {} K'.format(Tinf))
#py.xlabel('Time (s)')
#py.ylabel('Mass Fraction (dry basis)')
#py.legend(loc='best', numpoints=1)
#
#py.figure(3)
#py.plot(t, ligc, label='ligc')
#py.plot(t, g111, label='g1')
#py.plot(t, g222, label='g2')
#py.title('Lignin Reactions at T = {} K'.format(Tinf))
#py.xlabel('Time (s)')
#py.ylabel('Mass Fraction (dry basis)')
#py.legend(loc='best', numpoints=1)

py.show()
