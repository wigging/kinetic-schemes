# Contour plots of gas, tar, char production
# kinetic scheme based on Chan1985 and Blasi1993b

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
    
    return rww, rwg+rtg, rwt-rtg-rtc, rwc+rtc
    
    
#---- Parameters

rhow = 700  # wood density, kg/m^3    

tmax = 25   # max time, s
dt = 0.01   # time step, delta t
t = np.linspace(0, tmax, num=tmax/dt)   # time vector

#---- Evaluate kinetic rates

T = 673
rw673, rg673, rt673, rc673 = kinetics(rhow, T, tmax)

T = 723
rw723, rg723, rt723, rc723 = kinetics(rhow, T, tmax)

T = 773
rw773, rg773, rt773, rc773 = kinetics(rhow, T, tmax)

T = 823
rw823, rg823, rt823, rc823 = kinetics(rhow, T, tmax)

T = 873
rw873, rg873, rt873, rc873 = kinetics(rhow, T, tmax)

T = 923
rw923, rg923, rt923, rc923 = kinetics(rhow, T, tmax)

#---- Plot results

py.figure(1)
py.plot(t, rw773, label='wood')
py.plot(t, rg773, label='gas')
py.plot(t, rt773, label='tar')
py.plot(t, rc773, label='char')
py.xlabel('Time (s)')
py.ylabel('Rate ($kg/m^3s$)')
py.title(r'Kinetics at T = 773 K, $\rho_w$ = 700 $kg/m^3$')
py.legend(loc='best', numpoints=1)
py.grid()
py.show()

#---- Plot contour results

x = t
y = [673, 723, 773, 823]
z = np.array([rt673, rt723, rt773, rt823])

py.figure(2)
py.contourf(x, y, z, cmap='jet')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title(r'Tar Rate for $\rho_w$ = 700 $kg/m^3$')
py.colorbar().set_label('Tar Rate ($kg/m^3s$)')
py.grid()
py.show()

z = np.array([rg673, rg723, rg773, rg823])

py.figure(3)
py.contourf(x, y, z, cmap='jet')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title(r'Gas Rate for $\rho_w$ = 700 $kg/m^3$')
py.colorbar().set_label('Gas Rate ($kg/m^3s$)')
py.grid()
py.show()

z = np.array([rc673, rc723, rc773, rc823])

py.figure(4)
py.contourf(x, y, z, cmap='jet')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title(r'Char Rate for $\rho_w$ = 700 $kg/m^3$')
py.colorbar().set_label('Gas Rate ($kg/m^3s$)')
py.grid()
py.show()

