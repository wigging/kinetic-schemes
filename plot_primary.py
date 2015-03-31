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

# Parameters from Papadikis 2010a
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.01                               # time step, delta t
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

# convert concentrations to mass fraction for each scheme, (-)
wood0 = pw0/rhow;   gas0 = pg0/rhow;    tar0 = pt0/rhow;    char0 = pc0/rhow
wood1 = pw1/rhow;   gas1 = pg1/rhow;    tar1 = pt1/rhow;    char1 = pc1/rhow
wood2 = pw2/rhow;   gas2 = pg2/rhow;    tar2 = pt2/rhow;    char2 = pc2/rhow
wood3 = pw3/rhow;   gas3 = pg3/rhow;    tar3 = pt3/rhow;    char3 = pc0/rhow
wood4 = pw4/rhow;   gas4 = pg4/rhow;    tar4 = pt4/rhow;    char4 = pc4/rhow
wood5 = pw5/rhow;   gas5 = pg5/rhow;    tar5 = pt5/rhow;    char5 = pc5/rhow

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(1)
py.plot(t, wood0, label='Chan 1985')
py.plot(t, wood1, label='Font 1990 (fb)')
py.plot(t, wood2, label='Font 1990 (pp)')
py.plot(t, wood3, label='Janse 2000')
py.plot(t, wood4, label='Thurner 1981')
py.plot(t, wood5, label='Blasi 2001')
py.title('Wood Conversion, primary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Wood Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.figure(2)
py.plot(t, tar0, label='Chan 1985')
py.plot(t, tar1, label='Font 1990 (fb)')
py.plot(t, tar2, label='Font 1990 (pp)')
py.plot(t, tar3, label='Janse 2000')
py.plot(t, tar4, label='Thurner 1981')
py.plot(t, tar5, label='Blasi 2001')
py.title('Tar Yield, primary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Tar Mass Fraction (dry basis)')
py.legend(loc='best', numpoints=1)

py.show()
