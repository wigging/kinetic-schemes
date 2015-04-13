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
wtcell = 48
wthemi = 28
wtlig = 24

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.001                              # time step, delta t
tmax = 4                                # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Products from Cellulose, Hemicellulose, Lignin Reactions
#------------------------------------------------------------------------------

# arrays main groups and products as mass fractions, (-)
pmcell, pcell = kn.cell(rhow, wtcell, Tinf, dt, nt)     # cellulose
pmhemi, phemi = kn.hemi(rhow, wthemi, Tinf, dt, nt)     # hemicellulose
pmligc, pligc = kn.ligc(rhow, wtlig, Tinf, dt, nt)      # lignin-c
pmligh, pligh = kn.ligh(rhow, wtlig, Tinf, dt, nt)      # lignin-h
pmligo, pligo = kn.ligo(rhow, wtlig, Tinf, dt, nt)      # lignin-o

# main cellulose groups as mass fraction, (-)
cell = pmcell[0]
g1cell = pmcell[1]
cella = pmcell[2]
lvg = pmcell[3]
g4cell = pmcell[4]

# main hemicellulose groups as mass fraction, (-)
hemi = pmhemi[0]
g1hemi = pmhemi[1]
g2hemi = pmhemi[2]
g3hemi = pmhemi[3]
g4hemi = pmhemi[4]
xyl = pmhemi[5]

# main lignin-c groups as mass fraction, (-)
ligc = pmligc[0]
g1ligc = pmligc[1]
g2ligc = pmligc[2]

# main lignin-h groups as mass fraction, (-)
ligh = pmligh[0]
g1ligh = pmligh[1]
g2ligh = pmligh[2]
g3ligh = pmligh[3]
g4ligh = pmligh[4]
g5ligh = pmligh[5]
fe2macr1 = pmligh[6]

# main lignin-o groups as mass fraction, (-)
ligo = pmligo[0]
g1ligo = pmligo[1]
g2ligo = pmligo[2]
g3ligo = pmligo[3]
g4ligo = pmligo[4]
g5ligo = pmligo[5]
fe2macr2 = pmligo[6]

# main group mass balances as mass fraction, (-)
tcell = cell + g1cell + cella + lvg + g4cell            # cellulose
themi = hemi + g1hemi + g2hemi + g3hemi + g4hemi + xyl  # hemicellulose
tligc = ligc + g1ligc + g2ligc                                       # lignin-c
tligh = ligh + g1ligh + g2ligh + g3ligh + g4ligh + g5ligh + fe2macr1 # lignin-h
tligo = ligo + g1ligo + g2ligo + g3ligo + g4ligo + g5ligo + fe2macr2 # lignin-o

# Gas, Tar, Char from Cellulose, Hemicellulose, Lignin Reactions
#------------------------------------------------------------------------------

# chemical species as mass fraction, (-)
co = pcell[0] + phemi[0] + pligc[0] + pligh[0] + pligo[0]       # CO
co2 = pcell[1] + phemi[1] + pligc[1] + pligh[1] + pligo[1]      # CO2
ch2o = pcell[2] + phemi[2] + pligc[2] + pligh[2] + pligo[2]     # CH2O
hcooh = pcell[3] + phemi[3] + pligc[3] + pligh[3] + pligo[3]    # HCOOH
ch3oh = pcell[4] + phemi[4] + pligc[4] + pligh[4] + pligo[4]    # CH3OH
ch4 = pcell[5] + phemi[5] + pligc[5] + pligh[5] + pligo[5]      # CH4
glyox = pcell[6] + phemi[6] + pligc[6] + pligh[6] + pligo[6]    # Glyox (C2H2O2)
c2h4 = pcell[7] + phemi[7] + pligc[7] + pligh[7] + pligo[7]     # C2H4
c2h4o = pcell[8] + phemi[8] + pligc[8] + pligh[8] + pligo[8]    # C2H4O
haa = pcell[9] + phemi[9] + pligc[9] + pligh[9] + pligo[9]      # HAA (C2H4O2)
c2h5oh = pcell[10] + phemi[10] + pligc[10] + pligh[10] + pligo[10]  # C2H5OH
c3h6o = pcell[11] + phemi[11] + pligc[11] + pligh[11] + pligo[11]   # C3H6O
xyl = pcell[12] + phemi[12] + pligc[12] + pligh[12] + pligo[12]     # Xylose (C5H10O5)
c6h6o = pcell[13] + phemi[13] + pligc[13] + pligh[13] + pligo[13]   # C6H6O
hmfu = pcell[14] + phemi[14] + pligc[14] + pligh[14] + pligo[14]    # HMFU (C6H6O3)
lvg = pcell[15] + phemi[15] + pligc[15] + pligh[15] + pligo[15]     # LVG (C6H10O2)
coum = pcell[16] + phemi[16] + pligc[16] + pligh[16] + pligo[16]    # p-Coumaryl (C9H10O2)
fe2macr = pcell[17] + phemi[17] + pligc[17] + pligh[17] + pligo[17] # FE2MACR (C11H12O4)
h2 = pcell[18] + phemi[18] + pligc[18] + pligh[18] + pligo[18]      # H2
h2o = pcell[19] + phemi[19] + pligc[19] + pligh[19] + pligo[19]     # H2O
char = pcell[20] + phemi[20] + pligc[20] + pligh[20] + pligo[20]    # Char

# groups for gas and tar as mass fraction, (-)
gas = co + co2 + ch4 + c2h4 + h2
tar = ch2o + hcooh + ch3oh + glyox + c2h4o + haa + c2h5oh + c3h6o + xyl + \
      c6h6o + hmfu + lvg + coum + fe2macr

# chemical species mass balance as mass fraction, (-)
total = co + co2 + ch2o + hcooh + ch3oh + ch4 + glyox + c2h4 + c2h4o + haa + \
        c2h5oh + c3h6o + xyl + c6h6o + hmfu + lvg + coum + fe2macr + h2 + h2o + char


# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(1)
py.plot(t, cell, label='cell')
py.plot(t, cella, label='cella')
py.plot(t, g1cell, label='g1')
py.plot(t, g4cell, label='g4')
py.plot(t, lvg, label='lvg')
py.plot(t, tcell, '--k', label='total')
py.title('Cellulose Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(2)
py.plot(t, hemi, label='hemi')
py.plot(t, g1hemi, label='g1')
py.plot(t, g2hemi, label='g2')
py.plot(t, g3hemi, label='g3')
py.plot(t, g4hemi, label='g4')
py.plot(t, xyl, label='xyl')
py.plot(t, themi, '--k', label='total')
py.title('Hemicellulose Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(3)
py.plot(t, ligc, label='ligc')
py.plot(t, g1ligc, label='g1')
py.plot(t, g2ligc, label='g2')
py.plot(t, tligc, '--k', label='total')
py.title('Lignin-C Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(4)
py.plot(t, ligh, label='ligh')
py.plot(t, g1ligh, label='g1')
py.plot(t, g2ligh, label='g2')
py.plot(t, g3ligh, label='g3')
py.plot(t, g4ligh, label='g4')
py.plot(t, g5ligh, label='g5')
py.plot(t, fe2macr1, label='fe2macr')
py.plot(t, tligh, '--k', label='total')
py.title('Lignin-H Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(5)
py.plot(t, ligo, label='ligo')
py.plot(t, g1ligo, label='g1')
py.plot(t, g2ligo, label='g2')
py.plot(t, g3ligo, label='g3')
py.plot(t, g4ligo, label='g4')
py.plot(t, g5ligo, label='g5')
py.plot(t, fe2macr2, label='fe2macr')
py.plot(t, tligo, '--k', label='total')
py.title('Lignin-O Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(6)
py.plot(t, gas, label='gas')
py.plot(t, tar, label='tar')
py.plot(t, char, label='char')
py.plot(t, h2o, label='$\mathrm{H_{2}O}$')
py.plot(t, tar+h2o, label='tar + $\mathrm{H_{2}O}$')
py.title('Gas, Tar, Char at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Mass Fraction')
py.legend(loc='best', numpoints=1)

#py.figure(2)
#py.plot(t, h2o, label='$\mathrm{H_{2}O}$')
#py.plot(t, char, label='char')
#py.plot(t, haa, label='haa')
#py.plot(t, glyox, label='glyox')
#py.plot(t, c2h4o, label='$\mathrm{C_{2}H_{4}O}$')
#py.plot(t, hmfu, label='HMFU')
#py.plot(t, c3h6o, label='$\mathrm{C_{3}H_{6}O}$')
#py.plot(t, co2, ls='--', label='$\mathrm{CO_2}$')
#py.plot(t, h2, ls='--', label='$\mathrm{H_2}$')
#py.plot(t, ch2o, ls='--', label='$\mathrm{CH_2O}$')
#py.plot(t, co, ls='--', label='$\mathrm{CO}$')
#py.plot(t, ch4, ls='--', label='$\mathrm{CH_4}$')
#py.plot(t, hcooh, ls='--', label='$\mathrm{HCOOH}$')
#py.plot(t, lvg, ls='--', label='lvg')
#py.plot(t, total2, ls='--', label='total2')
#py.title('Cellulose Reactions at T = {} K'.format(Tinf))
#py.xlabel('Time (s)')
#py.ylabel('Mass Fraction (dry basis)')
#py.legend(loc='best', numpoints=1)

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
