"""
Plot product yields from primary & secondary reactions of biomass pyrolysis.
Reactions are provided by the kinetic scheme in a separate function file.

References:
See comments in the function file for references to a particular kinetic scheme.
"""

import numpy as np
import matplotlib.pyplot as py
import functions as fn

# Parameters
# ------------------------------------------------------------------------------

T = 773  # temperature for rate constants, K

dt = 0.005                              # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Products from Wood Kinetic Schemes
# ------------------------------------------------------------------------------

# store concentrations from primary reactions on a mass basis as kg/m^3
# row = concentration calculated from a particular kinetic scheme
# column = concentration at time step
wood = np.ones((13, nt))  # wood concentration array
tar = np.zeros((13, nt))  # tar concentration array

# products from primary reactions of Blasi 1993, Blasi 2001, Chan 1985,
# Font 1990, Janse 2000, Koufopanos 199, Liden 1988, Papadikis 2010,
# Sadhukhan 2009, Thurner 1981
for i in range(1, nt):
    wood[0, i], _, tar[0, i], _ = fn.blasi(wood[0, i-1], 0, tar[0, i-1], 0, T, dt, s=2)
    wood[1, i], _, tar[1, i], _ = fn.blasibranca(wood[1, i-1], 0, tar[1, i-1], 0, T, dt)
    wood[2, i], _, tar[2, i], _, _, _ = fn.chan(wood[2, i-1], 0, tar[2, i-1], 0, 0, 0, T, dt, s=2)
    wood[3, i], _, tar[3, i], _ = fn.font1(wood[3, i-1], 0, tar[3, i-1], 0, T, dt)
    wood[4, i], _, tar[4, i], _ = fn.font2(wood[4, i-1], 0, tar[4, i-1], 0, T, dt)
    wood[5, i], _, tar[5, i], _ = fn.janse(wood[5, i-1], 0, tar[5, i-1], 0, T, dt, s=2)
    wood[6, i], _, _, _, _ = fn.koufopanos(wood[6, i-1], 0, 0, 0, 0, T, dt, s=2)
    wood[7, i], _, tar[7, i], _ = fn.liden(wood[7, i-1], 0, tar[7, i-1], 0, T, dt, s=2)
    wood[8, i], _, tar[8, i], _ = fn.papadikis(wood[8, i-1], 0, tar[8, i-1], 0, T, dt, s=2)
    wood[9, i], _, _, _, _ = fn.sadhukhan(wood[9, i-1], 0, 0, 0, 0, T, dt, s=2)
    wood[10, i], _, tar[10, i], _ = fn.thurner(wood[10, i-1], 0, tar[10, i-1], 0, T, dt)

# Products from Ranzi 2014 Kinetic Scheme
# ------------------------------------------------------------------------------

# weight percent (%) cellulose, hemicellulose, lignin for beech wood
wtcell = 48
wthemi = 28
wtlig = 24

# arrays for Ranzi main groups and products as mass fractions, (-)
pmcell, pcell = fn.ranzicell(1, wtcell, T, dt, nt)    # cellulose
pmhemi, phemi = fn.ranzihemi(1, wthemi, T, dt, nt)    # hemicellulose
pmligc, pligc = fn.ranziligc(1, wtlig, T, dt, nt)     # lignin-c
pmligh, pligh = fn.ranziligh(1, wtlig, T, dt, nt)     # lignin-h
pmligo, pligo = fn.ranziligo(1, wtlig, T, dt, nt)     # lignin-o

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
wood_ranzi = pmcell[0] + pmhemi[0] + pmligc[0] + pmligh[0] + pmligo[0]
tar_ranzi = ch2o + hcooh + ch3oh + glyox + c2h4o + haa + c2h5oh + c3h6o + xyl + c6h6o + hmfu + lvg + coum + fe2macr

wood[11] = wood_ranzi
tar[11] = tar_ranzi

# Products from Miller and Bellan 1997 Kinetic Scheme
# ------------------------------------------------------------------------------

# composition of beech wood from Table 2 in paper
wtcell = 0.48   # cellulose mass fraction, (-)
wthemi = 0.28   # hemicellulose mass fraction, (-)
wtlig = 0.24    # lignin mass fraction, (-)

cella = np.ones(nt)*wtcell
hemia = np.ones(nt)*wthemi
liga = np.ones(nt)*wtlig

tar1, tar2, tar3 = np.zeros(nt), np.zeros(nt), np.zeros(nt)

for i in range(1, nt):
    cella[i], _, tar1[i], _ = fn.millercell_noR1(cella[i-1], 0, tar1[i-1], 0, T, dt, s=2)
    hemia[i], _, tar2[i], _ = fn.millerhemi_noR1(hemia[i-1], 0, tar2[i-1], 0, T, dt, s=2)
    liga[i], _, tar3[i], _ = fn.millerlig_noR1(liga[i-1], 0, tar3[i-1], 0, T, dt, s=2)

wood[12] = cella + hemia + liga
tar[12] = tar1 + tar2 + tar3

# Plot Results
# ------------------------------------------------------------------------------

py.ion()
py.close('all')

def despine():
    # remove top, right axis and tick marks
    ax = py.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    py.tick_params(axis='both', bottom='off', top='off', left='off', right='off')

py.figure(1)
py.plot(t, wood[7], lw=2, label='Liden 1988')
py.plot(t, wood[8], lw=2, label='Papadikis 2010')
py.plot(t, wood[5], lw=2, label='Janse 2000')
py.plot(t, wood[12], lw=2, label='Miller 1997')
py.title('Primary and secondary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Wood Concentration (mass fraction)')
py.legend(loc='best', numpoints=1, fontsize=11, frameon=False)
py.grid()
despine()

py.figure(2)
py.plot(t, tar[7], lw=2, label='Liden 1988')
py.plot(t, tar[8], lw=2, label='Papadikis 2010')
py.plot(t, tar[5], lw=2, label='Janse 2000')
py.plot(t, tar[12], lw=2, label='Miller 1997')
py.title('Primary and secondary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Tar Concentration (mass fraction)')
py.legend(loc='best', numpoints=1, fontsize=11, frameon=False)
py.grid()
despine()
