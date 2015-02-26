"""
Kinetic solution based on initial wood density and reaction rates, for example 
pw = pwi + rww*dt. Reactions are provided by a separate function file.

Function:
funcChanBlasi.py

Requirements:
Python 3, Numpy, Matplotlib

References:
1) Papadikis 2010a
2) Chan 1985 and Blasi 1993b
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
import funcChanBlasi as kn

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
pw = np.zeros(len(t))   # wood 
pg = np.zeros(len(t))   # gas
pt = np.zeros(len(t))   # tar
pc = np.zeros(len(t))   # char

pw[0] = rhow # initial wood density

# kinetics for primary and secondary reactions
for i in range(1, p):
    pw[i], pg[i], pt[i], pc[i] = kn.chanblasi(Tinf, pw, pg, pt, pc, dt, i)

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
 
py.figure(1)
py.plot(t, pw, label='wood')
py.plot(t, pg, label='gas')
py.plot(t, pt, label='tar')
py.plot(t, pc, label='char')
py.legend(loc='best', numpoints=1)
py.title('Reactions at T = %.f K' % Tinf)
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.grid()
py.show()
