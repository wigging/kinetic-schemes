"""
Batch reactor example with Debiagi 2018 biomass pyrolysis kinetics.

Results from this example should match those from the Cantera example. The
approach presented here uses SciPy `solve_ivp` to solve the system of ODEs
that represent the batch reactor material balance. However, it is tedious to
setup the function when dealing with a large system of reactions. A better
approach would be to use a stoichiometric matrix developed from a list of
reactions.

References
----------
P. Debiagi, G. Gentile, A. Cuoci, A. Frassoldati, E. Ranzi, and T. Faravelli.
A predictive model of biochar formation and characterization. Journal of
Analytical and Applied Pyrolysis, vol. 134, pp. 326-335, 2018.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
# -----------------------------------------------------------------------------

species = (
    'CELL', 'CELLA', 'CH2OHCHO', 'CHOCHO', 'CH3CHO', 'C6H6O3', 'C2H5CHO',
    'CH3OH', 'CH2O', 'CO', 'GCO', 'CO2', 'H2', 'H2O', 'GCH2O', 'HCOOH',
    'CH2OHCH2CHO', 'CH4', 'GH2', 'CHAR', 'C6H10O5', 'GCOH2', 'GMSW', 'HCE1',
    'HCE2'
)

# time vector to evaluate reaction rates [s]
time = np.linspace(0, 2.0)

# Batch reactor
# -----------------------------------------------------------------------------


def dcdt_debiagi(t, y):
    """
    Rate equations for Debiagi 2018 biomass pyrolysis kinetics.
    """
    R = 1.9859          # universal gas constant [cal/(mol K)]
    T = 773.15          # temperature [K]

    cell = y[0]         # cellulose mass fraction [-]
    cella = y[1]        # active cellulose mass fraction [-]
    gmsw = y[22]        # softwood hemicellulose mass fraction [-]

    # Cellulose reactions and rate constants
    # 1) CELL -> CELLA
    # 2) CELLA -> 0.40 CH2OHCHO + 0.03 CHOCHO + 0.17 CH3CHO + 0.25 C6H6O3 + 0.35 C2H5CHO + 0.20 CH3OH + 0.15 CH2O + 0.49 CO + 0.05 GCO + 0.43 CO2 + 0.13 H2 + 0.93 H2O + 0.05 GCH2O + 0.02 HCOOH + 0.05 CH2OHCH2CHO + 0.05 CH4 + 0.1 GH2 + 0.66 CHAR
    # 3) CELLA -> C6H10O5
    # 4) CELL -> 4.45 H2O + 5.45 CHAR + 0.12 GCOH2 + 0.18 GCH2O + 0.25 GCO + 0.125 GH2 + 0.125 H2
    k1 = 1.5e14 * np.exp(-47_000 / (R * T))
    k2 = 2.5e6 * np.exp(-19_100 / (R * T))
    k3 = 3.3 * T * np.exp(-10_000 / (R * T))
    k4 = 9e7 * np.exp(-31_000 / (R * T))

    # Hemicellulose reactions and rate constants
    # 5) GMSW -> 0.70 HCE1 + 0.30 HCE2
    k5 = 1e10 * np.exp(-31_000 / (R * T))

    # species reaction rate equations where r = dc/dt
    # mass fractions associated with each species are also given
    r_CELL = -k1 * cell
    r_CELLA = k1 * cell - k2 * cella - k3 * cella
    r_CH2OHCHO = k2 * cella * 0.1481
    r_CHOCHO = k2 * cella * 0.0107
    r_CH3CHO = k2 * cella * 0.0462
    r_C6H6O3 = k2 * cella * 0.1944
    r_C2H5CHO = k2 * cella * 0.1254
    r_CH3OH = k2 * cella * 0.0395
    r_CH2O = k2 * cella * 0.0278
    r_CO = k2 * cella * 0.0846
    r_GCO = k2 * cella * 0.008637 + k4 * cell * 0.04319
    r_CO2 = k2 * cella * 0.1167
    r_H2 = k2 * cella * 0.001616 + k4 * cell * 0.001554
    r_H2O = k2 * cella * 0.1033 + k4 * cell * 0.4944
    r_GCH2O = k2 * cella * 0.009259 + k4 * cell * 0.03333
    r_HCOOH = k2 * cella * 0.005677
    r_CH2OHCH2CHO = k2 * cella * 0.02284
    r_CH4 = k2 * cella * 0.004947
    r_GH2 = k2 * cella * 0.001243 + k4 * cell * 0.001554
    r_CHAR = k2 * cella * 0.04889 + k4 * cell * 0.4037
    r_C6H10O5 = k3 * cella
    r_GCOH2 = k4 * cell * 0.02222
    r_GMSW = -k5 * gmsw
    r_HCE1 = k5 * gmsw * 0.7
    r_HCE2 = k5 * gmsw * 0.3

    # system of ODEs
    dcdt = (
        r_CELL, r_CELLA, r_CH2OHCHO, r_CHOCHO, r_CH3CHO, r_C6H6O3, r_C2H5CHO,
        r_CH3OH, r_CH2O, r_CO, r_GCO, r_CO2, r_H2, r_H2O, r_GCH2O, r_HCOOH,
        r_CH2OHCH2CHO, r_CH4, r_GH2, r_CHAR, r_C6H10O5, r_GCOH2, r_GMSW, r_HCE1,
        r_HCE2
    )
    return dcdt


# initial mass fractions [-]
y0 = np.zeros(len(species))
y0[0] = 0.6
y0[22] = 0.4

# solve system of equations for a batch reactor
sol = solve_ivp(dcdt_debiagi, (time[0], time[-1]), y0, t_eval=time)

# Print
# ----------------------------------------------------------------------------

print('--- Initial ---')
print(f'CELL        {sol.y[0][0]}')
print(f'GMSW        {sol.y[22][0]}')

print('--- Final ---')
for i, sp in enumerate(species):
    print(f'{sp:11} {sol.y[i][-1]:.4f}')

# Plot
# -----------------------------------------------------------------------------


def config(ax, xlabel, ylabel):
    ax.grid(True, color='0.9')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax.set_frame_on(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(color='0.9')


fig, ax = plt.subplots(tight_layout=True)
for i, sp in enumerate(species):
    ax.plot(time, sol.y[i], label=sp)
config(ax, xlabel='Time [s]', ylabel='Mass fraction [-]')

plt.show()
