# Kinetic approach where wood density is calculated based on reaction rate
# pw = pwi + rww*dt
# kinetic scheme from Chan1985 and Blasi1993b

# use Python 3 print function
from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as py

py.close('all')

#--- parameters from Papadikis2010a
rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

#--- kinetic parameters from Chan1985 and Blasi1993b
# A = pre-exponential factor, 1/s and E = activation energy, kJ/mol
A1 = 1.3e8;  E1 = 140    # wood -> gas
A2 = 2e8;    E2 = 133    # wood -> tar
A3 = 1.08e7; E3 = 121    # wood -> char
A4 = 4.28e6; E4 = 108    # tar -> gas
A5 = 1e6;    E5 = 108    # tar -> char

R = 0.008314    # universal gas constant, kJ/mol*K

#--- kinetic equations
T = np.linspace(700, 1000)
K1 = A1 * np.exp(-E1 / (R * T))   # wood -> gas
K2 = A2 * np.exp(-E2 / (R * T))   # wood -> tar
K3 = A3 * np.exp(-E3 / (R * T))   # wood -> char
K4 = A4 * np.exp(-E4 / (R * T))   # tar -> gas
K5 = A5 * np.exp(-E5 / (R * T))   # tar -> char

#--- plot T vs K
py.figure(1)
py.plot(T, K1, label='K1 wood > gas')
py.plot(T, K2, label='K2 wood > tar')
py.plot(T, K3, label='K3 wood > char')
py.plot(T, K4, label='K4 tar > gas')
py.plot(T, K5, label='K5 tar > char')
py.legend(loc='best', numpoints=1)
py.xlabel('Temperature (K)')
py.ylabel('K, reaction rate constant (1/s)')
py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.grid()
py.show()

#--- initial calculations
dt = 0.01   # time step, delta t
tmax = 25   # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
p = len(t)  # total number of time steps

#--- kinetic equations
T = Tinf
K1 = A1 * np.exp(-E1 / (R * T))   # wood -> gas
K2 = A2 * np.exp(-E2 / (R * T))   # wood -> tar
K3 = A3 * np.exp(-E3 / (R * T))   # wood -> char
K4 = A4 * np.exp(-E4 / (R * T))   # tar -> gas
K5 = A5 * np.exp(-E5 / (R * T))   # tar -> char

#--- initial wood, gas, tar, char densities
pw = [rhow]
pg = [0]
pt = [0]
pc = [0]

rww = [-(K1 + K2 + K3) * pw[0]]  # rate of wood consumption
rwg = [K1 * pw[0]]   # rate of gas production from wood
rwt = [K2 * pw[0]]   # rate of tar production from wood
rwc = [K3 * pw[0]]   # rate of char production from wood
rtg = [K4 * pt[0]]   # rate of gas production from tar
rtc = [K5 * pt[0]]   # rate of char production from tar

#--- kinetics for primary and secondary reactions

for i in range(1, p):
    rww.append(-(K1 + K2 + K3) * pw[i-1])
    rwg.append(K1 * pw[i-1])
    rwt.append(K2 * pw[i-1])
    rwc.append(K3 * pw[i-1])
    rtg.append(K4 * pt[i-1])
    rtc.append(K5 * pt[i-1])
    pw.append(pw[i-1] + rww[i]*dt)
    pg.append(pg[i-1] + (rwg[i] + rtg[i])*dt)
    pt.append(pt[i-1] + (rwt[i] - rtg[i] - rtc[i])*dt)
    pc.append(pc[i-1] + (rwc[i] + rtc[i])*dt)

#--- plot t vs r and t vs p
py.figure(2)
py.plot(t, rww, label='rww')
py.plot(t, rwg, label='rwg')
py.plot(t, rwt, label='rwt')
py.plot(t, rwc, label='rwc')
py.plot(t, rtg, label='rtg')
py.plot(t, rtc, label='rtc')
py.legend(loc='best', numpoints=1)
py.xlabel('Time (s)')
py.ylabel('r, reaction rate (rho/s)')
py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.grid()
py.show()

py.figure(3)
py.plot(t, pw, label='pw')
py.plot(t, pg, label='pg')
py.plot(t, pt, label='pt')
py.plot(t, pc, label='pc')
py.legend(loc='best', numpoints=1)
py.xlabel('Time (s)')
py.ylabel('Density ($kg/m^3$)')
py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.grid()
py.show()
