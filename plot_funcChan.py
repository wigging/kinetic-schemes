"""
Plot Chan 1985 kinetic reactions.

Reference:
Chan, W.C.R., Kelbon, M. & Krieger, B.B., 1985. Modelling and experimental 
verification of physical and chemical processes during pyrolysis of a large 
biomass particle. Fuel, 64(11), pp.1505–1513.
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
import funcChan as kn

py.close('all')

# Parameters
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K
mc = 10     # moisture content, %

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.01   # time step, delta t
tmax = 25   # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Calculate Kinetic Reactions and Concentrations (Primary Reactions)
#------------------------------------------------------------------------------

# kinetics for primary reactions
pw, pg, pt, pc, pm, pv = kn.chan1(rhow, mc, Tinf, dt, nt)

# Calculate Kinetic Reactions and Concentrations (Primary + Secondary Reactions)
#------------------------------------------------------------------------------

# kinetics for primary and secondary reactions
pw2, pg2, pt2, pc2, pm2, pv2 = kn.chan3(rhow, mc, Tinf, dt, nt)

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
py.plot(t, pm, label='moisture')
py.plot(t, pv, label='water vapor')
py.title('Chan 1985 Primary Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.legend(loc='best', numpoints=1)
py.show()

py.figure(2)
py.plot(t, pw2, label='wood')
py.plot(t, pg2, label='gas')
py.plot(t, pt2, label='tar')
py.plot(t, pc2, label='char')
py.plot(t, pm2, label='moisture')
py.plot(t, pv2, label='water vapor')
py.title('Chan 1985 Primary & Secondary Reactions at T = {} K'.format(Tinf))
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.legend(loc='best', numpoints=1)
py.show()

