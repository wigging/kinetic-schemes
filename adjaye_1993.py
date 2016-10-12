"""
Plot products from Adjaye 1993 thesis for catalytic upgrading of bio-oil derived
from biomass pyrolysis. Kinetic parameters were taken from the thesis instead of
the 1995c and 1995d papers due to more realistic values given in the thesis.

References:
Adjaye, 1993. Thesis Appendix F, pp 311-319.
Adjaey, Bakhshi, 1995c. Biomass and Bioenergy, Vol. 8, No. 3, pp 131-149.
Adjaye, Bakhshi, 1995d. Biomass and Bioenergy, Vol. 8, No. 4, pp 265-277.
"""

import numpy as np
import matplotlib.pyplot as py

# Parameters
# -----------------------------------------------------------------------------

dt = 0.01                             # time step, delta t
tmax = 25                             # max time, s
t = np.linspace(0, tmax, num=tmax/dt) # time vector
nt = len(t)                           # total number of time steps


def adjaye(oil, nv, v, dt):
    # values for T = 410 degC
    Knv = 0.3
    Kr1 = 9.2e7
    Kcr = 4.1e-5
    Kc1 = 3.7e5
    Kv = 0.8
    Kd = 8.0e-4
    Ka = 6.1e-6
    Kg = 1.8e-4
    Kr2 = 37.0e5

    rOil = -(Knv + Kv)*oil
    rNV = Knv*oil - Kr1*nv - Kcr*nv - Kc1*nv
    rV = Kv*oil + Kcr*nv - Kd*v - Ka*v - Kg*v - Kr2*v
    nOil = oil + rOil*dt
    nNV = nv + rNV*dt
    nV = v + rV*dt

    return nOil, nNV, nV


oil = np.ones(nt)     # bio-oil concentration
nonvol = np.zeros(nt) # non-volatiles concentration
vol = np.zeros(nt)    # volatiles concentration

for i in range(1, nt):
    oil[i], nonvol[i], vol[i] = adjaye(oil[i-1], nonvol[i-1], vol[i-1], dt)

# Plot Results
# -----------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, oil, lw=2, label='oil')
py.plot(t, nonvol, lw=2, label='nonvol')
py.plot(t, vol, lw=2, label='vol')
py.title('Adjaye 1993 catalyst kinetics')
py.xlabel('Time')
py.ylabel('Concentration')
py.legend(loc='best', numpoints=1)
py.grid()
