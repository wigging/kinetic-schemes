# Kinetic approach where wood density is calculated based on reaction rate
# pw = pwi + rww*dt
# considers moisture content of the wood
# kinetic reactions scheme from Chan1985 and Blasi1993b

# use Python 3 print function
from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as py

# Parameters
#------------------------------------------------------------------------------

rhow = 700      # density of wood, kg/m^3
Tinf = 773      # ambient temp, K
mc = 0.20       # moisture content of wood
wf = 1-mc       # wood fraction in particle
R = 0.008314    # universal gas constant, kJ/mol*K

# Kinetic parameters from Chan1985 and Blasi1993b
# A = pre-exponential factor, 1/s and E = activation energy, kJ/mol
#-----------------------------------------------------------------------------

A1 = 1.3e8;     E1 = 140    # wood -> gas
A2 = 2e8;       E2 = 133    # wood -> tar
A3 = 1.08e7;    E3 = 121    # wood -> char
A4 = 4.28e6;    E4 = 108    # tar -> gas
A5 = 1e6;       E5 = 108    # tar -> char
Aw = 5.13e6;    Ew = 87.9   # water -> vapor

# Reaction rate constants to plot as function of temperature
#------------------------------------------------------------------------------

Tp = np.linspace(700, 1000)
K1p = A1 * np.exp(-E1 / (R * Tp))   # wood -> gas
K2p = A2 * np.exp(-E2 / (R * Tp))   # wood -> tar
K3p = A3 * np.exp(-E3 / (R * Tp))   # wood -> char
K4p = A4 * np.exp(-E4 / (R * Tp))   # tar -> gas
K5p = A5 * np.exp(-E5 / (R * Tp))   # tar -> char
Kwp = Aw * np.exp(-Ew / (R * Tp))   # water -> vapor

# Initial calculations
#------------------------------------------------------------------------------

dt = 0.01                               # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
p = len(t)                              # total number of time steps

# Reaction rate constants, 1/s
#------------------------------------------------------------------------------

T = Tinf                            # temperature
K1 = A1 * np.exp(-E1 / (R * T))     # wood -> gas
K2 = A2 * np.exp(-E2 / (R * T))     # wood -> tar
K3 = A3 * np.exp(-E3 / (R * T))     # wood -> char
K4 = A4 * np.exp(-E4 / (R * T))     # tar -> gas
K5 = A5 * np.exp(-E5 / (R * T))     # tar -> char
Kw = Aw * np.exp(-Ew / (R * T))     # water -> vapor

# Initial concentrations as density, p (rho) as  kg/m^3
#------------------------------------------------------------------------------

pw = [rhow*wf]  # wood conentration
pg = [0]        # gas concentration
pt = [0]        # tar concentration
pc = [0]        # char concentration
pwa = [rhow*mc] # water concentration
pva = [0]       # vapor concentration

# Initial reaction rates, r (rate) as rho/s
#------------------------------------------------------------------------------

rww = [-(K1 + K2 + K3) * pw[0]] # rate of wood consumption
rwg = [K1 * pw[0]]              # rate of gas production from wood
rwt = [K2 * pw[0]]              # rate of tar production from wood
rwc = [K3 * pw[0]]              # rate of char production from wood
rtg = [K4 * pt[0]]              # rate of gas production from tar
rtc = [K5 * pt[0]]              # rate of char production from tar
rwa = [-Kw * pwa[0]]            # rate of water vaporization
rva = [Kw * pwa[0]]             # rate of vapor production

# Calculate reaction rates & concentrations with primary & secondary reactions
#------------------------------------------------------------------------------

for i in range(1, p):
    rww.append(-(K1 + K2 + K3) * pw[i-1])
    rwg.append(K1 * pw[i-1])
    rwt.append(K2 * pw[i-1])
    rwc.append(K3 * pw[i-1])
    rtg.append(K4 * pt[i-1])
    rtc.append(K5 * pt[i-1])
    rwa.append(-Kw * pwa[i-1])
    rva.append(Kw * pwa[i-1])
    pw.append(pw[i-1] + rww[i]*dt)
    pg.append(pg[i-1] + (rwg[i] + rtg[i])*dt)
    pt.append(pt[i-1] + (rwt[i] - rtg[i] - rtc[i])*dt)
    pc.append(pc[i-1] + (rwc[i] + rtc[i])*dt)
    pwa.append(pwa[i-1] + rwa[i]*dt)
    pva.append(pva[i-1] + rva[i]*dt)


# Plot Results
#------------------------------------------------------------------------------

# close all previous plots
py.close('all')

# temperature vs reaction rate constant
py.figure(1)
py.plot(Tp, K1p, label=r'wood$\rightarrow$gas')
py.plot(Tp, K2p, label=r'wood$\rightarrow$tar')
py.plot(Tp, K3p, label=r'wood$\rightarrow$char')
py.plot(Tp, K4p, label=r'tar$\rightarrow$gas')
py.plot(Tp, K5p, label=r'tar$\rightarrow$char')
py.plot(Tp, Kwp, label=r'water$\rightarrow$vapor')
py.legend(loc='best', numpoints=1, fontsize=12)
py.xlabel('Temperature (K)')
py.ylabel('Reaction Rate Constant, K (1/s)')
py.title('Reaction Rate Constants')
py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.setp(py.gca().lines, lw=2)
py.grid()
py.show()

# time vs reaction rate
# rwa and rva are not displayed
py.figure(2)
py.plot(t, rww, label='wood')
py.plot(t, rwg, label=r'wood$\rightarrow$gas')
py.plot(t, rwt, label=r'wood$\rightarrow$tar')
py.plot(t, rwc, label=r'wood$\rightarrow$char')
py.plot(t, rtg, label=r'tar$\rightarrow$gas')
py.plot(t, rtc, label=r'tar$\rightarrow$char')
py.legend(loc='best', numpoints=1, fontsize=12)
py.xlabel('Time (s)')
py.ylabel('Reaction Rate, r (rho/s)')
py.title('Reaction Rates')
py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.setp(py.gca().lines, lw=2)
py.grid()
py.show()

# time vs density
py.figure(3)
py.plot(t, pw, 'green', label='wood')
py.plot(t, pg, 'orange', label='gas')
py.plot(t, pt, 'brown', label='tar')
py.plot(t, pc, 'black', label='char')
py.plot(t, pwa, 'blue', label='water')
py.plot(t, pva, 'cyan', label='vapor')
py.legend(loc='best', numpoints=1, fontsize=12)
py.xlabel('Time (s)')
py.ylabel('Density ($kg/m^3$)')
py.title('Wood Devol. and Product Yields')
py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.setp(py.gca().lines, lw=2)
py.grid()
py.show()
