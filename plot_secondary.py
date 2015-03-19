"""
Plot product yields from primary & secondary reactions of biomass pyrolysis. 
Reactions are provided by the kinetic scheme in a separate function file.

Function Files:
funcPapadikis.py
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

import funcPapadikis as kn0
import funcJanse as kn1
import funcMiller as kn2

# Parameters from Papadikis 2010a
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.0001                               # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
p = len(t)                              # total number of time steps

# Calculate Kinetic Reactions and Concentrations - Groups
#------------------------------------------------------------------------------

# arrays for wood, gas, tar, char concentrations as a density, kg/m^3
# row = concentration for a particular kinetic scheme
# column = time step
pw = np.zeros([3, p])   # wood array
pg = np.zeros([3, p])   # gas array
pt = np.zeros([3, p])   # tar array
pc = np.zeros([3, p])   # char array

pw[:] = rhow # initial wood density

# kinetics for primary & secondary reactions - groups
pw[0], pg[0], pt[0], pc[0] = kn0.papadikis2(Tinf, pw[0], dt, p)
pw[1], pg[1], pt[1], pc[1] = kn1.janse2(Tinf, pw[1], dt, p)

# Calculate Kinetic Reactions and Concentrations - Cell, Hemi, Lig
#------------------------------------------------------------------------------

# array of product concentrations as a density, kg/m^3
pcell, pcella, pgas1, ptar1, pchar1 = kn2.cell(Tinf, pw[2], dt, p)
phemi, phemia, pgas2, ptar2, pchar2 = kn2.hemi(Tinf, pw[2], dt, p)
plig, pliga, pgas3, ptar3, pchar3 = kn2.lig(Tinf, pw[2], dt, p)

# sum components to wood, gas, tar, char groups, kg/m^3
pw[2] = pcell + phemi + plig
pactive = pcella + phemia + pliga
pg[2] = pgas1 + pgas2 + pgas3
pt[2] = ptar1 + ptar2 + ptar3
pc[2] = pchar1 + pchar2 + pchar3

# convert concentrations to percent
wood = pw/rhow*100
gas = pg/rhow*100
tar = pt/rhow*100
char = pc/rhow*100

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(3)
py.plot(t, wood[0], label='Papadikis 2010')
py.plot(t, wood[1], label='Janse 2000')
py.plot(t, wood[2], label='Miller 1997')
py.title('Wood Conversion, primary & secondary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Wood Conversion (% Dry Basis)')
py.legend(loc='best', numpoints=1)
py.show()

py.figure(4)
py.plot(t, tar[0], label='Papadikis 2010')
py.plot(t, tar[1], label='Janse 2000')
py.plot(t, tar[2], label='Miller 1997')
py.title('Tar Yield, primary & secondary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Tar Yield (% Dry Basis)')
py.legend(loc='best', numpoints=1)
py.show()
