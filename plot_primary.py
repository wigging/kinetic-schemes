"""
Plot product yields from primary reactions of biomass pyrolysis. Reactions are 
provided by the kinetic scheme in a separate function file.

Function Files:
funcChan.py
funcFont.py
funcThurner.py
funcJanse.py

Requirements:
Python 3, Numpy, Matplotlib

References:
See comments in each function file.
"""

# Modules and Function Files
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
py.close('all')

import funcChan as kn0
import funcFont as kn1
import funcJanse as kn2
import funcThurner as kn3
import funcBlasiBranca as kn4
import funcRanzi as kn5

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
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Calculate Kinetic Reactions and Concentrations
#------------------------------------------------------------------------------

# kinetics for primary reactions for product concentrations, kg/m^3
pw0, pg0, pt0, pc0 = kn0.chan2(rhow, Tinf, dt, nt)
pw1, pg1, pt1, pc1 = kn1.font1(rhow, Tinf, dt, nt)
pw2, pg2, pt2, pc2 = kn1.font2(rhow, Tinf, dt, nt)
pw3, pg3, pt3, pc3 = kn2.janse1(rhow, Tinf, dt, nt)
pw4, pg4, pt4, pc4 = kn3.thurner(rhow, Tinf, dt, nt)
pw5, pg5, pt5, pc5 = kn4.blasibranca(rhow, Tinf, dt, nt)

# assign concentrations as mass fraction for each scheme, (-)
wood0 = pw0/rhow;   gas0 = pg0/rhow;    tar0 = pt0/rhow;    char0 = pc0/rhow
wood1 = pw1/rhow;   gas1 = pg1/rhow;    tar1 = pt1/rhow;    char1 = pc1/rhow
wood2 = pw2/rhow;   gas2 = pg2/rhow;    tar2 = pt2/rhow;    char2 = pc2/rhow
wood3 = pw3/rhow;   gas3 = pg3/rhow;    tar3 = pt3/rhow;    char3 = pc0/rhow
wood4 = pw4/rhow;   gas4 = pg4/rhow;    tar4 = pt4/rhow;    char4 = pc4/rhow
wood5 = pw5/rhow;   gas5 = pg5/rhow;    tar5 = pt5/rhow;    char5 = pc5/rhow

# Calculate Kinetic Reactions and Concentrations from Ranzi
#------------------------------------------------------------------------------

# arrays for Ranzi main groups and products as mass fractions, (-)
pmcell, pcell = kn5.cell(rhow, wtcell, Tinf, dt, nt)     # cellulose
pmhemi, phemi = kn5.hemi(rhow, wthemi, Tinf, dt, nt)     # hemicellulose
pmligc, pligc = kn5.ligc(rhow, wtlig, Tinf, dt, nt)      # lignin-c
pmligh, pligh = kn5.ligh(rhow, wtlig, Tinf, dt, nt)      # lignin-h
pmligo, pligo = kn5.ligo(rhow, wtlig, Tinf, dt, nt)      # lignin-o

# chemical species from Ranzi as mass fraction, (-)
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

# groups from Ranzi for wood and tar as mass fraction, (-)
wood6 = pmcell[0] + pmhemi[0] + pmligc[0] + pmligh[0] + pmligo[0]
tar = ch2o + hcooh + ch3oh + glyox + c2h4o + haa + c2h5oh + c3h6o + xyl + \
      c6h6o + hmfu + lvg + coum + fe2macr

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['xtick.major.size'] = 0
py.rcParams['ytick.major.pad'] = 8
py.rcParams['ytick.major.size'] = 0
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True
py.rcParams['legend.framealpha'] = 0

py.figure(1)
py.plot(t, wood0, label='Chan 1985')
py.plot(t, wood1, label='Font 1990 (fb)')
py.plot(t, wood2, label='Font 1990 (pp)')
py.plot(t, wood3, label='Janse 2000')
py.plot(t, wood4, label='Thurner 1981')
py.plot(t, wood5, label='Blasi 2001')
py.plot(t, wood6, label='Ranzi 2014')
py.title('Wood Conversion, primary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Wood Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)
py.gca().spines['right'].set_visible(False)
py.gca().spines['top'].set_visible(False)

py.figure(2)
py.plot(t, tar0, label='Chan 1985')
py.plot(t, tar1, label='Font 1990 (fb)')
py.plot(t, tar2, label='Font 1990 (pp)')
py.plot(t, tar3, label='Janse 2000')
py.plot(t, tar4, label='Thurner 1981')
py.plot(t, tar5, label='Blasi 2001')
py.plot(t, tar, '--k', label='Ranzi 2014 (tar)')
py.plot(t, tar+h2o, '-k', label='Ranzi 2014 (tar+$\mathrm{H_{2}O}$)')
py.title('Tar Yield, primary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Tar Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)
py.gca().spines['right'].set_visible(False)
py.gca().spines['top'].set_visible(False)

py.show()
