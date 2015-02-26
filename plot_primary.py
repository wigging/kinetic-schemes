"""
Plot product yields from primary reactions of biomass pyrolysis. Reactions are 
provided by the kinetic scheme in a separate function file.

Function Files:
funcChanBlasi.py
funcThurnerMann.py
funcBlasiBranca.py
funcFont.py

Requirements:
Python 3, Numpy, Matplotlib

References:
See comments function file.
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
import funcChanBlasi as kn1
import funcThurnerMann as kn2
import funcBlasiBranca as kn3
import funcFont as kn4

py.close('all')

# Parameters from Papadikis 2010a
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.01   # time step, delta t
tmax = 25   # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
p = len(t)  # total number of time steps

# Calculate Kinetic Reactions and Concentrations
#------------------------------------------------------------------------------

# vectors for wood, char, gas, tar concentrations as a density, kg/m^3
pw = np.zeros([5, len(t)])   # wood 
pc = np.zeros([5, len(t)])   # char
pg = np.zeros([5, len(t)])   # gas
pt = np.zeros([5, len(t)])   # tar 

pw[:] = rhow # initial wood density

# kinetics for primary and secondary reactions
pw1 = pw[0]; pc1 = pc[0]; pg1 = pg[0]; pt1 = pt[0]
pw2 = pw[1]; pc2 = pc[1]; pg2 = pg[1]; pt2 = pt[1]
pw3 = pw[2]; pc3 = pc[2]; pg3 = pg[2]; pt3 = pt[2]
pw4 = pw[3]; pc4 = pc[3]; pg4 = pg[3]; pt4 = pt[3]
pw5 = pw[4]; pc5 = pc[4]; pg5 = pg[4]; pt5 = pt[4]

# kinetics for primary and secondary reactions
for i in range(1, p):
    pw[0,i], pg[0,i], pt[0,i], pc[0,i] = kn1.chanblasi_p(Tinf, pw1, pg1, pt1, pc1, dt, i)
    pw[1,i], pg[1,i], pt[1,i], pc[1,i] = kn2.thurnermann(Tinf, pw2, pg2, pt2, pc2, dt, i)
    pw[2,i], pg[2,i], pt[2,i], pc[2,i] = kn3.blasibranca(Tinf, pw3, pg3, pt3, pc3, dt, i)
    pw[3,i], pg[3,i], pt[3,i], pc[3,i] = kn4.font1(Tinf, pw4, pg4, pt4, pc4, dt, i)
    pw[4,i], pg[4,i], pt[4,i], pc[4,i] = kn4.font2(Tinf, pw5, pg5, pt5, pc5, dt, i)

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(1)
py.plot(t, pw[0], label='ChanBlasi')
py.plot(t, pw[1], label='ThurnerMann')
py.plot(t, pw[2], label='BlasiBranca')
py.plot(t, pw[3], label='Font1')
py.plot(t, pw[4], label='Font2')
py.title('Wood Conversion at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.legend(loc='best', numpoints=1)
py.show()

py.figure(2)
py.plot(t, pt[0], label='ChanBlasi')
py.plot(t, pt[1], label='ThurnerMann')
py.plot(t, pt[2], label='BlasiBranca')
py.plot(t, pt[3], label='Font1')
py.plot(t, pt[4], label='Font2')
py.title('Tar Yield at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.legend(loc='best', numpoints=1)
py.show()
