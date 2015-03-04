"""
Analytical solution for a system of kinetic reactions for biomass pyrolysis.
Wood density as function of time, for example pw = rhow*exp(-(K1+K2+K3)*t).

Requirements:
Pyton 3, Numpy, Matplotlib

References:
1) Papadikis, Gu, Bridgwater, 2010. Fuel Processing Technology, 91(1), pp.68–79.
2) Chan, Kelbon, Krieger, 1985. Fuel, 64(11), pp.1505–1513.
3) Liden, Berruti, Scott, 1988. Chemical Engineering Communications, 65, pp.207–221.
4) Blasi, 1993. Combustion Science and Technology, 90, pp.315–340.
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py

# Parameters from Papadikis 2010a
#------------------------------------------------------------------------------
rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# Kinetic parameters from Chan 1985, Liden 1988, Blasi 1993b
#------------------------------------------------------------------------------

# A = pre-exponential factor, 1/s and E = activation energy, kJ/mol
A1 = 1.3e8;  E1 = 140    # wood -> gas
A2 = 2e8;    E2 = 133    # wood -> tar
A3 = 1.08e7; E3 = 121    # wood -> char
A4 = 4.28e6; E4 = 108    # tar -> gas
A5 = 1e6;    E5 = 108    # tar -> char

R = 0.008314    # universal gas constant, kJ/mol*K

# Initial calculations
#------------------------------------------------------------------------------

dt = 0.01   # time step, delta t
tmax = 25   # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
p = len(t)  # total number of time steps

# Calculate Kinetic Reactions and Concentrations
#------------------------------------------------------------------------------

T = Tinf
K1 = A1 * np.exp(-E1 / (R * T))   # wood -> gas
K2 = A2 * np.exp(-E2 / (R * T))   # wood -> tar
K3 = A3 * np.exp(-E3 / (R * T))   # wood -> char
K4 = A4 * np.exp(-E4 / (R * T))   # tar -> gas
K5 = A5 * np.exp(-E5 / (R * T))   # tar -> char

# wood concentration as a density vs time, kg/m^3
pw = rhow * np.exp(-(K1 + K2 + K3) * t)

# initial gas, tar, char concentrations as a density, kg/m^3
pg = [0]    # gas
pt = [0]    # tar
pc = [0]    # char

rwg = [K1 * pw[0]]   # rate of gas production from wood
rwt = [K2 * pw[0]]   # rate of tar production from wood
rwc = [K3 * pw[0]]   # rate of char production from wood
rtg = [K4 * pt[0]]   # rate of gas production from tar
rtc = [K5 * pt[0]]   # rate of char production from tar

# kinetics for primary and secondary reactions
for i in range(1, p):
    rwg.append(K1 * pw[i-1])
    rwt.append(K2 * pw[i-1])
    rwc.append(K3 * pw[i-1])
    rtg.append(K4 * pt[i-1])
    rtc.append(K5 * pt[i-1])
    pg.append(pg[i-1] + (rwg[i] + rtg[i])*dt)
    pt.append(pt[i-1] + (rwt[i] - rtg[i] - rtc[i])*dt)
    pc.append(pc[i-1] + (rwc[i] + rtc[i])*dt)

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
py.legend(loc='best', numpoints=1)
py.title('Analytical Solution, reactions at T = %.f K' % Tinf)
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.show()
