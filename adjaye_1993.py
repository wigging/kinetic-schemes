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

    # concentrations
    oil = c[0]      # bio-oil concentration
    nonvol = c[1]   # non-volatiles concentration
    vol = c[2]      # volatiles concentration
    # coke = c[3]     # coke concentration
    # res = c[4]      # residuce concentration
    # aqueous = c[4]  # aqueous concentration
    # organic = c[6]  # organic distillate concentration
    # hydro = c[7]    # hydrocarbon concentration
    # gas = c[8]      # gas concentration

    # rate constants
    Knv = 0.4       # bio-oil -> non-volatiles
    Kvo = 1.1       # bio-oil -> volatiles
    Kcr = 6.8e-5    # non-volatiles -> volatiles
    Kc1 = 3.4e-5    # non-volatiles -> coke
    Kr1 = 3.3e-7    # non-volatiles -> residue
    Kr2 = 8.3e-5    # volatiles -> residue
    Kaq = 3.4e-6    # volatiles -> aqueous fraction
    Kod = 8.3e-4    # volatiles -> organic distillates
    Kga = 1.1e-4    # volatiles -> gas
    Kc2 = 6.4e-5    # volatiles -> coke
    Khc = 6.0e-4    # volatiles -> hydrocarbons

    # reaction orders
    nv = 1      # bio-oil -> non-volatiles
    vo = 1      # bio-oil -> volatiles
    cr = 0.9    # non-volatiles -> volatiles
    c1 = 0.9    # non-volatiles -> coke
    r1 = 2.2    # non-volatiles -> residue
    r2 = 1      # volatiles -> residue
    aq = 1.6    # volatiles -> aqueous fraction
    od = 1      # volatiles -> organic distillates
    ga = 0.7    # volatiles -> gas
    c2 = 1.1    # volatiles -> coke
    hc = 0.9    # volatiles -> hydrocarbons

    # reaction rates
    rOil = -Knv*oil**nv - Kvo*oil**vo
    rNonVol = Knv*oil**nv - Kcr*nonvol**cr - Kc1*nonvol**c1 - Kr1*nonvol**r1
    rVol = Kvo*oil**vo + Kcr*nonvol**cr - Kr2*vol**r2 - Kc2*vol**c2 - Kaq*vol**aq - Kod*vol**od - Khc*vol**hc - Kga*vol**ga
    rCoke = Kc1*nonvol**c1 + Kc2*vol**c2
    rRes = Kr1*nonvol**r1 + Kr2*vol**r2
    rAq = Kaq*vol**aq
    rOrg = Kod*vol**od
    rHydro = Khc*vol**hc
    rGas = Kga*vol**ga

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
plt.xlabel('Time')
plt.ylabel('Concentration')
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
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(loc='best', numpoints=1)
plt.grid()

