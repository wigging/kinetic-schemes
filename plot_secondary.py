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

# Parameters from Papadikis 2010a
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.01                               # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
p = len(t)                              # total number of time steps

# Calculate Kinetic Reactions and Concentrations
#------------------------------------------------------------------------------

# arrays for wood, gas, tar, char concentrations as a density, kg/m^3
# row = concentration for a particular kinetic scheme
# column = time step
pw = np.zeros([5, len(t)])   # wood array
pg = np.zeros([5, len(t)])   # gas array
pt = np.zeros([5, len(t)])   # tar array
pc = np.zeros([5, len(t)])   # char array

pw[:] = rhow # initial wood density

# assign kinetic scheme to a particular row
pw0 = pw[0]; pc0 = pc[0]; pg0 = pg[0]; pt0 = pt[0]  # Papadikis kinetic scheme
pw1 = pw[1]; pc1 = pc[1]; pg1 = pg[1]; pt1 = pt[1]  # Janse kinetic scheme


# assign kinetic scheme to a particular row

# kinetics for primary reactions
for i in range(1, p):
    pw[0,i], pg[0,i], pt[0,i], pc[0,i] = kn0.papadikis2(Tinf, pw0, pg0, pt0, pc0, dt, i)
    pw[1,i], pg[1,i], pt[1,i], pc[1,i] = kn1.janse2(Tinf, pw1, pg1, pt1, pc1, dt, i)

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
py.title('Wood Conversion, primary & secondary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Wood Conversion (% Dry Basis)')
py.legend(loc='best', numpoints=1)
py.show()

py.figure(4)
py.plot(t, tar[0], label='Papadikis 2010')
py.plot(t, tar[1], label='Janse 2000')
py.title('Tar Yield, primary & secondary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Tar Yield (% Dry Basis)')
py.legend(loc='best', numpoints=1)
py.show()
