"""
Adjaye kinetics with adjusted kinetic parameters from Vivek at NREL.

References:

Adjaye, J.D. & Bakhshi, N.N., 1995. Catalytic Conversion of a Biomass-Derived
Oil to Fuels and Chemicals I: Model Compound Studies and Reaction Pathways.
Biomass & bioenergy, 8(3), pp.131–149.

Adjaye, J.D. & Bakhshi, N.N., 1995. Catalytic Conversion of a Biomass-Derived
Oil to Fuels and Chemicals II: Chemical Kinetics, Parameter Estimation and
Model Predictions. Biomass and Bioenergy, 8(4), pp.265–277.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp

# Parameters
# -----------------------------------------------------------------------------

dt = 0.01                               # time step, delta t
tmax = 1000                             # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Approach 1 - Euler method
# -----------------------------------------------------------------------------

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

# Approach 2 - SciPy odeint solver with primary and secondary reactions
# -----------------------------------------------------------------------------

def adjaye2(c, t):
    """
    Adjaye kinetics at 370C using Scipy odeint solver.

    Inputs
    ------
    c = concentration vector where c[0] = oil, c[1] = nonvol, c[2] = vol
    t = time vector

    Returns
    -------
    r = rates vector
    """

    # concentrations as kmol/m^3
    Coil = c[0]     # bio-oil
    Cnon = c[1]     # non-volatiles
    Cvol = c[2]     # volatiles

    # rate constants as K, m^3/(kg hr)
    # reaction orders as n, (-)
    K1 = 0.4        # bio-oil -> non-volatiles
    n1 = 1.0
    K2 = 1.1        # bio-oil -> volatiles
    n2 = 1
    K3 = 6.8e-5     # non-volatiles -> volatiles
    n3 = 0.9
    K4 = 3.4e-5     # non-volatiles -> coke
    n4 = 0.9
    K5 = 3.3e-7     # non-volatiles -> residue
    n5 = 2.2
    K6 = 8.3e-5     # volatiles -> residue
    n6 = 1
    K7 = 3.4e-6     # volatiles -> aqueous fraction
    n7 = 1.6
    K8 = 8.3e-4     # volatiles -> organic distillates
    n8 = 1
    K9 = 1.1e-4     # volatiles -> gas
    n9 = 0.7
    K10 = 6.4e-5    # volatiles -> coke
    n10 = 1.1
    K11 = 6.0e-4    # volatiles -> hydrocarbons
    n11 = 0.9

    # reaction rates
    rOil = -K1*Coil**n1 - K2*Coil**n2
    rNonVol = K1*Coil**n1 - K3*Cnon**n3 - K4*Cnon**n4 - K5*Cnon**n5
    rVol = (K2*Coil**n2 + K3*Cnon**n3 - K6*Cvol**n6 - K7*Cvol**n7 - K8*Cvol**n8
            - K9*Cvol**n9 - K10*Cvol**n10 - K11*Cvol**n11)
    rCoke = K4*Cnon**n4 + K10*Cvol**n10
    rRes = K5*Cnon**n5 + K6*Cvol**n6
    rAq = K7*Cvol**n7
    rOrg = K8*Cvol**n8
    rHydro = K11*Cvol**n11
    rGas = K9*Cvol**n9

    # return list of reaction rates
    dcdt = [rOil, rNonVol, rVol, rCoke, rRes, rAq, rOrg, rHydro, rGas]

    return dcdt


conc = sp.odeint(adjaye2, [1, 0, 0, 0, 0, 0, 0, 0, 0], t)

# Plot
# ------------------------------------------------------------------------------

plt.ion()
plt.close('all')

plt.figure(1)
plt.plot(t, oil, lw=2, label='oil')
plt.plot(t, nonvol, lw=2, label='nonvol')
plt.plot(t, vol, lw=2, label='vol')
plt.xlabel('Time (hr)')
plt.ylabel('Concentration (kmol/m^3)')
plt.legend(loc='best', numpoints=1)
plt.grid()

plt.figure(2)
plt.plot(t, conc[:, 0], lw=2, label='oil')
plt.plot(t, conc[:, 1], lw=2, label='nonvol')
plt.plot(t, conc[:, 2], lw=2, label='vol')
plt.plot(t, conc[:, 3], lw=2, label='coke')
plt.plot(t, conc[:, 4], lw=2, label='res')
plt.plot(t, conc[:, 5], lw=2, label='aq')
plt.plot(t, conc[:, 6], lw=2, label='org')
plt.plot(t, conc[:, 7], lw=2, label='hydro')
plt.plot(t, conc[:, 8], lw=2, label='gas')
plt.xlabel('Time (hr)')
plt.ylabel('Concentration (kmol/m^3)')
plt.legend(loc='best', numpoints=1)
plt.grid()
# plt.savefig('/Users/Gavin/Desktop/adjaye.pdf', bbox_inches='tight')
