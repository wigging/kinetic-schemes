"""
Plot primary and secondary kinetic reactions for biomass pyrolysis.
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
import funcJanse as kn

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

# Calculate Kinetic Reactions and Concentrations (Primary Reactions)
#------------------------------------------------------------------------------

# vectors for wood, gas, tar, char concentrations as a density, kg/m^3
pw = np.zeros(len(t))   # wood 
pg = np.zeros(len(t))   # gas
pt = np.zeros(len(t))   # tar
pc = np.zeros(len(t))   # char

pw[:] = rhow        # initial wood density

# kinetics for primary reactions
for i in range(1, p):
    pw[i], pg[i], pt[i], pc[i] = kn.janse1(Tinf, pw, pg, pt, pc, dt, i)
    
# Calculate Kinetic Reactions and Concentrations (Primary + Secondary Reactions)
#------------------------------------------------------------------------------

# vectors for wood, gas, tar, char concentrations as a density, kg/m^3
pw2 = np.zeros(len(t))   # wood 
pg2 = np.zeros(len(t))   # gas
pt2 = np.zeros(len(t))   # tar
pc2 = np.zeros(len(t))   # char

pw2[:] = rhow        # initial wood density

# kinetics for primary and secondary reactions
for i in range(1, p):
    pw2[i], pg2[i], pt2[i], pc2[i] = kn.janse2(Tinf, pw2, pg2, pt2, pc2, dt, i)

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(1)
py.plot(t, pw, label='wood')
py.plot(t, pg, label='gas')
py.plot(t, pt, label='tar')
py.plot(t, pc, label='char')
py.title('Primary Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.legend(loc='best', numpoints=1)
py.show()

py.figure(2)
py.plot(t, pw2, label='wood')
py.plot(t, pg2, label='gas')
py.plot(t, pt2, label='tar')
py.plot(t, pc2, label='char')
py.title('Primary & Secondary Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.legend(loc='best', numpoints=1)
py.show()

