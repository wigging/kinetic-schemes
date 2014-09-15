# Contour Plots for kinetic scheme based on Chan1985 and Blasi1993b

# use Python 3 functions
from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as py

py.close('all')

#---- Kinetics function

def kinetics(rhow, T, tmax):
    A1 = 1.3e8;  E1 = 140    # wood -> gas
    A2 = 2e8;    E2 = 133    # wood -> tar
    A3 = 1.08e7; E3 = 121    # wood -> char
    A4 = 4.28e6; E4 = 108    # tar -> gas
    A5 = 1e6;    E5 = 108    # tar -> char
    R = 0.008314             # universal gas constant, kJ/mol*K
    
    dt = 0.01                               # time step, delta t
    t = np.linspace(0, tmax, num=tmax/dt)   # time vector
    p = len(t)                              # total number of time steps
    
    K1 = A1 * np.exp(-E1 / (R * T))   # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))   # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))   # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T))   # tar -> gas
    K5 = A5 * np.exp(-E5 / (R * T))   # tar -> char
    
    pw = np.zeros(p)
    pg = np.zeros(p)
    pt = np.zeros(p)
    pc = np.zeros(p)

    rww = np.zeros(p)   # rate of wood consumption
    rwg = np.zeros(p)   # rate of gas production from wood
    rwt = np.zeros(p)   # rate of tar production from wood
    rwc = np.zeros(p)   # rate of char production from wood
    rtg = np.zeros(p)   # rate of gas production from tar
    rtc = np.zeros(p)   # rate of char production from tar

    pw[0] = rhow
    rww[0] = -(K1 + K2 + K3) * pw[0]
    rwg[0] = K1 * pw[0]
    rwt[0] = K2 * pw[0]
    rwc[0] = K3 * pw[0]
    rtg[0] = K4 * pt[0]
    rtc[0] = K5 * pt[0]
    
    for i in range(1, p):
        rww[i] = -(K1 + K2 + K3) * pw[i-1]
        rwg[i] = K1 * pw[i-1]
        rwt[i] = K2 * pw[i-1]
        rwc[i] = K3 * pw[i-1]
        rtg[i] = K4 * pt[i-1]
        rtc[i] = K5 * pt[i-1]
        pw[i] = pw[i-1] + rww[i]*dt
        pg[i] = pg[i-1] + (rwg[i] + rtg[i])*dt
        pt[i] = pt[i-1] + (rwt[i] - rtg[i] - rtc[i])*dt
        pc[i] = pc[i-1] + (rwc[i] + rtc[i])*dt
    
    return pw, pg, pt, pc

#---- Parameters

rhow = 700  # wood density, kg/m^3    

tmax = 25   # max time, s
dt = 0.01   # time step, delta t
t = np.linspace(0, tmax, num=tmax/dt)   # time vector

#---- Evaluate kinetics at different temperatures

T = 673
pw673, pg673, pt673, pc673 = kinetics(rhow, T, tmax)

T = 723
pw723, pg723, pt723, pc723 = kinetics(rhow, T, tmax)

T = 773
pw773, pg773, pt773, pc773 = kinetics(rhow, T, tmax)

T = 823
pw823, pg823, pt823, pc823 = kinetics(rhow, T, tmax)

T = 873
pw873, pg873, pt873, pc873 = kinetics(rhow, T, tmax)

T = 923
pw923, pg923, pt923, pc923 = kinetics(rhow, T, tmax)

#---- Plot as density

x = t
y = [673, 723, 773, 823, 873, 923]

z = np.array([pt673, pt723, pt773, pt823, pt873, pt923])

py.figure(1)
py.contourf(x, y, z, cmap='jet')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title(r'Tar Production, $\rho_w$ = 700 $kg/m^3$')
py.colorbar().set_label('Tar ($kg/m^3$)')
py.show()

#---- Plot as normalized to rhow

z = np.array([pt673/rhow, pt723/rhow, pt773/rhow, pt823/rhow, pt873/rhow, pt923/rhow])

py.figure(2)
py.contourf(x, y, z, cmap='jet')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title(r'Tar Production, $\rho_w$ = 700 $kg/m^3$')
py.colorbar().set_label('Tar (-)')
py.show()

z = np.array([pg673/rhow, pg723/rhow, pg773/rhow, pg823/rhow, pg873/rhow, pg923/rhow])

py.figure(3)
py.contourf(x, y, z, cmap='jet')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title(r'Gas Production, $\rho_w$ = 700 $kg/m^3$')
py.colorbar().set_label('Gas (-)')
py.show()

z = np.array([pc673/rhow, pc723/rhow, pc773/rhow, pc823/rhow, pc873/rhow, pc923/rhow])

py.figure(4)
py.contourf(x, y, z, cmap='jet')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title(r'Char Production, $\rho_w$ = 700 $kg/m^3$')
py.colorbar().set_label('Char (-)')
py.show()

z = np.array([pw673/rhow, pw723/rhow, pw773/rhow, pw823/rhow, pw873/rhow, pw923/rhow])

py.figure(5)
py.contourf(x, y, z, cmap='jet')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title(r'Wood Consumption, $\rho_w$ = 700 $kg/m^3$')
py.colorbar().set_label('Wood (-)')
py.show()

#---- Plot all as normalized to rhow

py.figure(6)

py.subplot(311)
z = np.array([pt673/rhow, pt723/rhow, pt773/rhow, pt823/rhow, pt873/rhow, pt923/rhow])
py.contourf(x, y, z, cmap='jet')
py.colorbar().set_label('Tar (-)')
py.ylabel('Temperature (K)')
py.title(r'Chan1985 & Blasi1993b Kinetics, $\rho_w$ = 700 $kg/m^3$')

py.subplot(312)
z = np.array([pg673/rhow, pg723/rhow, pg773/rhow, pg823/rhow, pg873/rhow, pg923/rhow])
py.contourf(x, y, z, cmap='jet')
py.colorbar().set_label('Gas (-)')
py.ylabel('Temperature (K)')

py.subplot(313)
z = np.array([pc673/rhow, pc723/rhow, pc773/rhow, pc823/rhow, pc873/rhow, pc923/rhow])
py.contourf(x, y, z, cmap='jet')
py.colorbar().set_label('Char (-)')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')

py.show()
