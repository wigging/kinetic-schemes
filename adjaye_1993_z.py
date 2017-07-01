"""
Adjaye kinetics adjusted using parameters from Vivek
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
# -----------------------------------------------------------------------------

dt = 0.01                               # time step, delta t
tmax = 85                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

def adjaye(oil, nonvol, vol, dt):
    """adjaye kinetics for 370 degC"""
    Knv = 0.4
    Kv = 1.1
    b1 = 1
    b2 = 1

    Kcr = 6.8e-5
    cr = 0.9
    Kr1 = 3.3e-7
    r1 = 2.2
    Kc1 = 3.4e-5
    c1 = 0.9

    Kr2 = 8.3e-5
    r2 = 1
    Kc2 = 6.4e-5
    c2 = 1.1
    Ka = 3.4e-6
    a = 1.6
    Kd = 8.3e-4
    d = 1
    Kh = 6e-4
    h = 0.9
    Kg = 1.1e-4
    g = 0.7

    rOil = -Knv*oil**b1 - Kv*oil**b2
    nOil = oil + rOil*dt

    rNonVol = Knv*oil**b1 - Kcr*nonvol**cr - Kr1*nonvol**r1 - Kc1*nonvol**c1
    nNonVol = nonvol + rNonVol*dt

    rVol = Kv*oil**b2 + Kcr*nonvol**cr - Kr2*vol**r2 - Kc2*vol**c2 - Ka*vol**a - Kd*vol**d - Kh*vol**h - Kg*vol**g
    nVol = vol + rVol*dt

    return nOil, nNonVol, nVol


oil = np.ones(nt)       # bio-oil concentration
nonvol = np.zeros(nt)   # non-volatiles concentration
vol = np.zeros(nt)      # volatiles concentration

for i in range(1, nt):
    oil[i], nonvol[i], vol[i] = adjaye(oil[i-1], nonvol[i-1], vol[i-1], dt)

# Plot
# ------------------------------------------------------------------------------

plt.ion()
plt.close('all')

plt.figure(1)
plt.plot(t, oil, lw=2, label='oil')
plt.plot(t, nonvol, lw=2, label='nonvol')
plt.plot(t, vol, lw=2, label='vol')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='best', numpoints=1)
plt.grid()

