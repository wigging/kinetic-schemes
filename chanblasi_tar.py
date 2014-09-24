# Tar production at different wood densities
# kinetic scheme based on Chan1985 and Blasi1993b

# use Python 3 functions
from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as py

py.close('all')

# Function for kinetic scheme
#-----------------------------------------------------------------------------

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
    

# Effect of density on tar production at T = 723 K
#-----------------------------------------------------------------------------
    
# parameters

T = 723     # temperature, K

tmax = 25   # max time, s
dt = 0.01   # time step, delta t
t = np.linspace(0, tmax, num=tmax/dt)   # time vector

# evaluate at wood density rhow = 400 kg/m^3

rhow = 400
pw400, pg400, pt400, pc400 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 500 kg/m^3

rhow = 500
pw500, pg500, pt500, pc500 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 600 kg/m^3

rhow = 600
pw600, pg600, pt600, pc600 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 700 kg/m^3

rhow = 700
pw700, pg700, pt700, pc700 = kinetics(rhow, T, tmax)

# plot results

py.figure(1)
py.plot(pw400, pt400, label='400')
py.plot(pw500, pt500, label='500')
py.plot(pw600, pt600, label='600')
py.plot(pw700, pt700, label='700')
py.xlabel('Wood ($kg/m^3$)')
py.ylabel('Tar ($kg/m^3$)')
py.ylim([0, 250])
py.title('Tar vs Wood Conversion at T = %.fK' % T)
py.legend(loc='best', numpoints=1)
py.gca().invert_xaxis()
py.grid()
py.show()

# Effect of density on tar production at T = 773 K
#-----------------------------------------------------------------------------
    
# parameters

T = 773     # temperature, K

tmax = 25   # max time, s
dt = 0.01   # time step, delta t
t = np.linspace(0, tmax, num=tmax/dt)   # time vector

# evaluate at wood density rhow = 400 kg/m^3

rhow = 400
pw400, pg400, pt400, pc400 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 500 kg/m^3

rhow = 500
pw500, pg500, pt500, pc500 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 600 kg/m^3

rhow = 600
pw600, pg600, pt600, pc600 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 700 kg/m^3

rhow = 700
pw700, pg700, pt700, pc700 = kinetics(rhow, T, tmax)

# plot results

py.figure(2)
py.plot(pw400, pt400, label='400')
py.plot(pw500, pt500, label='500')
py.plot(pw600, pt600, label='600')
py.plot(pw700, pt700, label='700')
py.xlabel('Wood ($kg/m^3$)')
py.ylabel('Tar ($kg/m^3$)')
py.ylim([0, 250])
py.title('Tar vs Wood Conversion at T = %.fK' % T)
py.legend(loc='best', numpoints=1)
py.gca().invert_xaxis()
py.grid()
py.show()

# Effect of density on tar production at T = 823 K
#-----------------------------------------------------------------------------
    
# parameters

T = 823     # temperature, K

tmax = 25   # max time, s
dt = 0.01   # time step, delta t
t = np.linspace(0, tmax, num=tmax/dt)   # time vector

# evaluate at wood density rhow = 400 kg/m^3

rhow = 400
pw400, pg400, pt400, pc400 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 500 kg/m^3

rhow = 500
pw500, pg500, pt500, pc500 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 600 kg/m^3

rhow = 600
pw600, pg600, pt600, pc600 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 700 kg/m^3

rhow = 700
pw700, pg700, pt700, pc700 = kinetics(rhow, T, tmax)

# plot results

py.figure(3)
py.plot(pw400, pt400, label='400')
py.plot(pw500, pt500, label='500')
py.plot(pw600, pt600, label='600')
py.plot(pw700, pt700, label='700')
py.xlabel('Wood ($kg/m^3$)')
py.ylabel('Tar ($kg/m^3$)')
py.title('Tar vs Wood Conversion at T = %.fK' % T)
py.legend(loc='best', numpoints=1)
py.gca().invert_xaxis()
py.grid()
py.show()

# Effect of density on tar production at T = 773 K and t = 25 s
#-----------------------------------------------------------------------------

# parameters

T = 773     # temperature, K

tmax = 25   # max time, s
dt = 0.01   # time step, delta t
t = np.linspace(0, tmax, num=tmax/dt)   # time vector

# evaluate at wood density rhow = 300 kg/m^3
rhow = 300
pw300, pg300, pt300, pc300 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 400 kg/m^3
rhow = 400
pw400, pg400, pt400, pc400 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 500 kg/m^3
rhow = 500
pw500, pg500, pt500, pc500 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 600 kg/m^3
rhow = 600
pw600, pg600, pt600, pc600 = kinetics(rhow, T, tmax)

# evaluate at wood density rhow = 700 kg/m^3
rhow = 700
pw700, pg700, pt700, pc700 = kinetics(rhow, T, tmax)

# plot results

py.figure(4)
py.plot(t, pt300, label='300 $kg/m^3$')
py.plot(t, pt400, label='400 $kg/m^3$')
py.plot(t, pt500, label='500 $kg/m^3$')
py.plot(t, pt600, label='600 $kg/m^3$')
py.plot(t, pt700, label='700 $kg/m^3$')
py.fill_between(t, 0, pt300, facecolor= 'b')
py.fill_between(t, pt300, pt400, facecolor='g')
py.fill_between(t, pt400, pt500, facecolor='r')
py.fill_between(t, pt500, pt600, facecolor='c')
py.fill_between(t, pt600, pt700, facecolor='m')
py.xlabel('Time (s)')
py.ylabel('Tar ($kg/m^3$)')
py.ylim([0, 200])
py.title('Tar vs Time at T = 773 K at t = 25 s')
py.legend(loc='best', numpoints=1)
py.grid()
py.show()

# Effect of temperature on tar production for rho = 700 kg/m^3 and t = 25 s
#-----------------------------------------------------------------------------

# parameters
rhow = 700  # wood density, kg/m^3    

tmax = 25   # max time, s
dt = 0.01   # time step, delta t
t = np.linspace(0, tmax, num=tmax/dt)   # time vector

# evaluate at temperature = 773 K
T = 773
pw773, pg773, pt773, pc773 = kinetics(rhow, T, tmax)

# evaluate at temperature = 823 K
T = 823
pw823, pg823, pt823, pc823 = kinetics(rhow, T, tmax)

# plot results
py.figure(5)
py.plot(pw773, pt773, label='773 K')
py.plot(pw823, pt823, label='823 K')
py.xlabel('Wood ($kg/m^3$)')
py.ylabel('Tar ($kg/m^3$)')
py.title('Tar vs Wood Conversion for rho = 700 $kg/m^3$ and t = 25 s')
py.legend(loc='best', numpoints=1)
py.gca().invert_xaxis()
py.grid()
py.show()

