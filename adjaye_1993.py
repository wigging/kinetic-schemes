"""
Adjaye kinetic scheme for conversion of bio-oil over zeolite catalyst. Reaction
pathways were adjusted based on available kinetic parameters. See function
definition for more details.

References:

Adjaye, J.D. Catalytic Conversion of Biomass-Derived Oils to Fuels and
Chemicals. Thesis submitted to Department of Chemical Engineering, University
of Saskatchewan, 1993.

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

t1 = np.linspace(0, 10, num=1000)     # time vector for 10 seconds
t2 = np.linspace(0, 1800, num=1000)   # time vector for 1800 seconds or 30 minutes

# Concentrations from Adjaye
# -----------------------------------------------------------------------------

def adjaye(c, t):
    """
    Calculate concentrations based on Adjaye kinetic scheme for bio-oil
    conversion in fixed bed of zeolite catalyst. Kinetic parameters were
    evaluated at 370C. Reaction pathways adjusted based on available data in
    Adjaye 1995 papers.

    Inputs
    ------
    c = concentration vector, kmol/m^3
        where c[0] = oil, c[1] = nonvol, c[2] = vol
    t = time vector, s

    Returns
    -------
    r = reaction rates vector, kmol/(kg cat. hr)
    """

    # concentrations as kmol/m^3
    Coil = c[0]     # bio-oil
    Cnon = c[1]     # non-volatiles
    Cvol = c[2]     # volatiles

    # rate constants as K, assumed to be m^3/(kg hr)
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

    # reaction rates as kmol/(kg cat. hr)
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


co = [1, 0, 0, 0, 0, 0, 0, 0, 0]    # initial concentrations
conc1 = sp.odeint(adjaye, co, t1)   # concentrations for times t1
conc2 = sp.odeint(adjaye, co, t2)   # concentration for times t2

# check mass balance, values in array should be same as initial concentration
total1 = conc1.sum(axis=1)
total2 = conc2.sum(axis=1)

# Plot
# ------------------------------------------------------------------------------

plt.ion()
plt.close('all')

plt.figure(1)
plt.plot(t1, conc1[:, 0], lw=2, label='oil')
plt.plot(t1, conc1[:, 1], lw=2, label='nonvol')
plt.plot(t1, conc1[:, 2], lw=2, label='vol')
plt.plot(t1, conc1[:, 3], lw=2, label='coke')
plt.plot(t1, conc1[:, 4], lw=2, label='res')
plt.plot(t1, conc1[:, 5], lw=2, label='aq')
plt.plot(t1, conc1[:, 6], lw=2, label='org')
plt.plot(t1, conc1[:, 7], lw=2, label='hydro')
plt.plot(t1, conc1[:, 8], lw=2, label='gas')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (kmol/m^3)')
plt.legend(fontsize=10, frameon=False, loc='best', numpoints=1)
plt.grid()
# plt.savefig('/Users/Gavin/Desktop/adjaye1.pdf', bbox_inches='tight')

plt.figure(2)
plt.plot(t2, conc2[:, 0], lw=2, label='oil')
plt.plot(t2, conc2[:, 1], lw=2, label='nonvol')
plt.plot(t2, conc2[:, 2], lw=2, label='vol')
plt.plot(t2, conc2[:, 3], lw=2, label='coke')
plt.plot(t2, conc2[:, 4], lw=2, label='res')
plt.plot(t2, conc2[:, 5], lw=2, label='aq')
plt.plot(t2, conc2[:, 6], lw=2, label='org')
plt.plot(t2, conc2[:, 7], lw=2, label='hydro')
plt.plot(t2, conc2[:, 8], lw=2, label='gas')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (kmol/m^3)')
plt.legend(fontsize=10, frameon=False, loc='best', numpoints=1)
plt.grid()
# plt.ylim(ymin=0)
# plt.savefig('/Users/Gavin/Desktop/adjaye2.pdf', bbox_inches='tight')

