"""
Plot product yields from primary & secondary reactions of biomass pyrolysis. 
Reactions are provided by the kinetic scheme in a separate function file.

Function Files:
funcPapadikis.py
funcJanse.py
funcMiller.py

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

import funcPapadikis as kn0
import funcJanse as kn1
import funcMiller as kn2

# Parameters from Papadikis 2010a
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.0001                             # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Calculate Kinetic Reactions and Concentrations - Groups
#------------------------------------------------------------------------------

# kinetics for product concentrations as a density, kg/m^3
pw0, pg0, pt0, pc0 = kn0.papadikis2(rhow, Tinf, dt, nt)
pw1, pg1, pt1, pc1 = kn1.janse2(rhow, Tinf, dt, nt)

# Calculate Kinetic Reactions and Concentrations - Cell, Hemi, Lig
#------------------------------------------------------------------------------

# kinetics for product concentrations as a density, kg/m^3
pcell, pcella, pgas1, ptar1, pchar1 = kn2.cell(rhow, Tinf, dt, nt)
phemi, phemia, pgas2, ptar2, pchar2 = kn2.hemi(rhow, Tinf, dt, nt)
plig, pliga, pgas3, ptar3, pchar3 = kn2.lig(rhow, Tinf, dt, nt)

# sum components to wood, gas, tar, char groups, kg/m^3
pactive = pcella + phemia + pliga
pw2 = pcell + phemi + plig + pactive
pg2 = pgas1 + pgas2 + pgas3
pt2 = ptar1 + ptar2 + ptar3
pc2 = pchar1 + pchar2 + pchar3

# Convert concentrations to mass fraction for each scheme
#------------------------------------------------------------------------------

wood0 = pw0/rhow;   gas0 = pg0/rhow;    tar0 = pt0/rhow;    char0 = pc0/rhow
wood1 = pw1/rhow;   gas1 = pg1/rhow;    tar1 = pt1/rhow;    char1 = pc1/rhow
wood2 = pw2/rhow;   gas2 = pg2/rhow;    tar2 = pt2/rhow;    char2 = pc2/rhow

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(3)
py.plot(t, wood0, label='Papadikis 2010')
py.plot(t, wood1, label='Janse 2000')
py.plot(t, wood2, label='Miller 1997')
py.title('Wood Conversion, primary & secondary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Wood Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(4)
py.plot(t, tar0, label='Papadikis 2010')
py.plot(t, tar1, label='Janse 2000')
py.plot(t, tar2, label='Miller 1997')
py.title('Tar Yield, primary & secondary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Tar Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.show()
