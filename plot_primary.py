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
pw0 = pw[0]; pc0 = pc[0]; pg0 = pg[0]; pt0 = pt[0]  # Chan kinetic scheme
pw1 = pw[1]; pc1 = pc[1]; pg1 = pg[1]; pt1 = pt[1]  # Font scheme (fluidized bed)
pw2 = pw[2]; pc2 = pc[2]; pg2 = pg[2]; pt2 = pt[2]  # Font scheme (pyroprobe)
pw3 = pw[3]; pc3 = pc[3]; pg3 = pg[3]; pt3 = pt[3]  # Janse kinetic scheme
pw4 = pw[4]; pc4 = pc[4]; pg4 = pg[4]; pt4 = pt[4]  # Thurner kinetic scheme

# assign kinetic scheme to a particular row

# kinetics for primary reactions
for i in range(1, p):
    pw[0,i], pg[0,i], pt[0,i], pc[0,i] = kn0.chan2(Tinf, pw0, pg0, pt0, pc0, dt, i)
    pw[1,i], pg[1,i], pt[1,i], pc[1,i] = kn1.font1(Tinf, pw1, pg1, pt1, pc1, dt, i)
    pw[2,i], pg[2,i], pt[2,i], pc[2,i] = kn1.font2(Tinf, pw2, pg2, pt2, pc2, dt, i)
    pw[3,i], pg[3,i], pt[3,i], pc[3,i] = kn2.janse1(Tinf, pw3, pg3, pt3, pc3, dt, i)
    pw[4,i], pg[4,i], pt[4,i], pc[4,i] = kn3.thurner(Tinf, pw4, pg4, pt4, pc4, dt, i)

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

py.figure(1)
py.plot(t, wood[0], label='Chan 1985')
py.plot(t, wood[1], label='Font 1990 (fb)')
py.plot(t, wood[2], label='Font 1990 (pp)')
py.plot(t, wood[3], label='Janse 2000')
py.plot(t, wood[4], label='Thurner 1981')
py.title('Wood Conversion, primary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Wood Conversion (% Dry Basis)')
py.legend(loc='best', numpoints=1)
py.show()

py.figure(2)
py.plot(t, tar[0], label='Chan 1985')
py.plot(t, tar[1], label='Font 1990 (fb)')
py.plot(t, tar[2], label='Font 1990 (pp)')
py.plot(t, tar[3], label='Janse 2000')
py.plot(t, tar[4], label='Thurner 1981')
py.title('Tar Yield, primary reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Tar Yield (% Dry Basis)')
py.legend(loc='best', numpoints=1)
py.show()
