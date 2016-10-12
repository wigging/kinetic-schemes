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
import scipy.integrate as sp
import matplotlib.pyplot as py

# Parameters
# -----------------------------------------------------------------------------

dt = 0.0001                           # time step, delta t
tmax = 25                             # max time, s
t = np.linspace(0, tmax, num=tmax/dt) # time vector
nt = len(t)                           # total number of time steps

# Approach 1 - Euler method with primary reactions only
# -----------------------------------------------------------------------------

def adjaye(oil, nv, v, dt):
    # values for T = 410 degC
    Knv = 0.3
    Kv = 0.8
    rOil = -(Knv + Kv)*oil
    rNV = Knv*oil
    rV = Kv*oil
    nOil = oil + rOil*dt
    nNV = nv + rNV*dt
    nV = v + rV*dt
    return nOil, nNV, nV

oil = np.ones(nt)     # bio-oil concentration
nonvol = np.zeros(nt) # non-volatiles concentration
vol = np.zeros(nt)    # volatiles concentration

for i in range(1, nt):
    oil[i], nonvol[i], vol[i] = adjaye(oil[i-1], nonvol[i-1], vol[i-1], dt)

# Approach 2 - SciPy odeint solver with primary reactions only
# -----------------------------------------------------------------------------

def rates(c, t):
    """
    w = wood-oil as conc[0]
    nv = non-volatiles as conc[1]
    v = volatiles as conc[2]
    """
    Knv = 0.3
    Kv = 0.8
    rw = -(Knv + Kv) * c[0]
    rnv = Knv * c[0]
    rv = Kv * c[0]
    return [rw, rnv, rv]

cc = sp.odeint(rates, [1, 0, 0], t)
# cc[cc < 0] = 0  # make negative values zero

# Approach 3 - SciPy odeint solver with primary and secondary reactions
# -----------------------------------------------------------------------------

def dCdt(c, t):
    """
    w = wood-oil as conc[0]
    nv = non-volatiles as conc[1]
    v = volatiles as conc[2]
    """
    Knv = 0.3
    Kv = 0.8
    Kr1 = 9.2e7;    r1 = 2.5
    Kcr = 4.1e-5;   cr = 0.9
    Kc1 = 3.7e5;    c1 = 1.1
    Kd = 8.0e-4;    d = 0.9
    Ka = 6.1e-6;    a = 1.4
    Kg = 1.8e-4;    g = 0.8
    Kr2 = 37.0e5;   r2 = 0.7
    rw = -(Knv + Kv)*c[0]
    rnv = Knv*c[0] - Kr1*c[1]**r1 - Kcr*c[1]**cr - Kc1*c[1]**c1
    rv = Kv*c[0] + Kcr*c[1]**cr - Kd*c[2]**d - Ka*c[2]**a - Kg*c[2]**g - Kr2*c[2]**r2
    return [rw, rnv, rv]

cc2 = sp.odeint(dCdt, [1, 0, 0], t)

# Approach 4 - SciPy runge kutta solver with primary and secondary reactions
# -----------------------------------------------------------------------------

def dcdt(t, c):
    """
    w = wood-oil as conc[0]
    nv = non-volatiles as conc[1]
    v = volatiles as conc[2]
    """
    Knv = 0.3
    Kv = 0.8
    Kr1 = 9.2e7;    r1 = 2.5
    Kcr = 4.1e-5;   cr = 0.9
    Kc1 = 3.7e5;    c1 = 1.1
    Kd = 8.0e-4;    d = 0.9
    Ka = 6.1e-6;    a = 1.4
    Kg = 1.8e-4;    g = 0.8
    Kr2 = 37.0e5;   r2 = 0.7
    rw = -(Knv + Kv)*c[0]
    rnv = Knv*c[0] - Kr1*c[1]**r1 - Kcr*c[1]**cr - Kc1*c[1]**c1
    rv = Kv*c[0] + Kcr*c[1]**cr - Kd*c[2]**d - Ka*c[2]**a - Kg*c[2]**g - Kr2*c[2]**r2
    return [rw, rnv, rv]

# store concentrations
Coil = np.ones(nt)      # bio-oil concentration
Cnvol = np.zeros(nt)    # non-volatiles concentration
Cvol = np.zeros(nt)     # volatiles concentration

# Setup the ode integrator where 'dopri5' is Runge-Kutta 4th order with 'bdf'
# for backward differentiation formula
r = sp.ode(dcdt).set_integrator('dopri5', nsteps=10000)
r.set_initial_value([1, 0, 0], 0)

# integrate the odes for each time step then store the results
k = 1
while r.successful() and r.t < tmax-dt:
    r.integrate(r.t+dt)
    Coil[k] = r.y[0]
    Cnvol[k] = r.y[1]
    Cvol[k] = r.y[2]
    k += 1

# Plot Results
# -----------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, oil, lw=2, label='oil')
py.plot(t, nonvol, lw=2, label='nonvol')
py.plot(t, vol, lw=2, label='vol')
py.title('Euler method, primary reactions')
py.xlabel('Time')
py.ylabel('Concentration')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t, cc[:, 0], lw=2, label='wood-oil')
py.plot(t, cc[:, 1], lw=2, label='non-vol')
py.plot(t, cc[:, 2], lw=2, label='vol')
py.title('SciPy odeint, primary reactions')
py.xlabel('Time')
py.ylabel('Concentration')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(3)
py.plot(t, cc2[:, 0], lw=2, label='wood-oil')
py.plot(t, cc2[:, 1], lw=2, label='non-vol')
py.plot(t, cc2[:, 2], lw=2, label='vol')
py.title('SciPy odeint, primary + secondary reactions')
py.xlabel('Time')
py.ylabel('Concentration')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(4)
py.plot(t, Coil, lw=2, label='oil')
py.plot(t, Cnvol, lw=2, label='nonvol')
py.plot(t, Cvol, lw=2, label='vol')
py.title('SciPy ode, primary + secondary reactions')
py.xlabel('Time')
py.ylabel('Concentration')
py.legend(loc='best', numpoints=1)
py.grid()

