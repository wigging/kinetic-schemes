"""
Plot primary and secondary kinetic reactions of biomass pyrolysis for the 
functional group (wood, gas, tar, char) kinetic schemes.
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
#import funcJanse as kn
import funcPapadikis as kn

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
nt = len(t)  # total number of time steps

# Calculate Kinetic Reactions and Concentrations (Primary Reactions)
#------------------------------------------------------------------------------

# kinetics for primary reactions    
pw1, pg1, pt1, pc1 = kn.papadikis1(rhow, Tinf, dt, nt)
    
# Calculate Kinetic Reactions and Concentrations (Primary + Secondary Reactions)
#------------------------------------------------------------------------------

# kinetics for primary and secondary reactions
pw2, pg2, pt2, pc2 = kn.papadikis2(rhow, Tinf, dt, nt)

# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(1)
py.plot(t, pw1, label='wood')
py.plot(t, pg1, label='gas')
py.plot(t, pt1, label='tar')
py.plot(t, pc1, label='char')
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

