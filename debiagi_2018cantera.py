"""
Cantera batch reactor example with Debiagi 2018 biomass pyrolysis kinetics.
Requires the `debiagi_sw.cti` file which defines the species and kinetic
reactions associated with softwood pyrolysis.

References
----------
P. Debiagi, G. Gentile, A. Cuoci, A. Frassoldati, E. Ranzi, and T. Faravelli.
A predictive model of biochar formation and characterization. Journal of
Analytical and Applied Pyrolysis, vol. 134, pp. 326-335, 2018.
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# Parameters
# ----------------------------------------------------------------------------

tk = 773.15     # reactor temperature [K]
p = 101325.0    # reactor pressure [Pa]

# initial mass fractions [-]
y0 = 'CELL:0.482 GMSW:0.262 LIGO:0.001 LIGH:0.173 LIGC:0.026 TANN:0.000 TGL:0.049'

# time vector to evaluate reaction rates [s]
time = np.linspace(0, 2.0, 100)

# Cantera batch reactor
# ----------------------------------------------------------------------------

gas = ct.Solution('debiagi_sw.cti')
gas.TPY = tk, p, y0
r = ct.IdealGasReactor(gas, energy='off')

sim = ct.ReactorNet([r])
states = ct.SolutionArray(gas, extra=['t'])

for tm in time:
    sim.advance(tm)
    states.append(r.thermo.state, t=tm)

# Print
# ----------------------------------------------------------------------------

print(f"""
--- Initial ---
CELL        {states('CELL').Y[:, 0][0]:.4f}
GMSW        {states('GMSW').Y[:, 0][0]:.4f}
LIGO        {states('LIGO').Y[:, 0][0]:.4f}
LIGH        {states('LIGH').Y[:, 0][0]:.4f}
LIGC        {states('LIGC').Y[:, 0][0]:.4f}
TANN        {states('TANN').Y[:, 0][0]:.4f}
TGL         {states('TGL').Y[:, 0][0]:.4f}
""")

print('--- Final ---')
for sp in states.species_names:
    print(f"{sp:11} {states(sp).Y[:, 0][-1]:.4f}")

# Plot
# ----------------------------------------------------------------------------


def config(ax, xlabel, ylabel):
    ax.grid(True, color='0.9')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    ax.set_frame_on(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(color='0.9')


fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states('CELL').Y[:, 0], label='CELL')
ax.plot(states.t, states('GMSW').Y[:, 0], label='GMSW')
ax.plot(states.t, states('LIGC').Y[:, 0], label='LIGC')
ax.plot(states.t, states('LIGH').Y[:, 0], label='LIGH')
ax.plot(states.t, states('LIGO').Y[:, 0], label='LIGO')
ax.plot(states.t, states('TANN').Y[:, 0], label='TANN')
ax.plot(states.t, states('TGL').Y[:, 0], label='TGL')
config(ax, xlabel='Time [s]', ylabel='Mass fraction [-]')

fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states('CELLA').Y[:, 0], label='CELLA')
ax.plot(states.t, states('HCE1').Y[:, 0], label='HCE1')
ax.plot(states.t, states('HCE2').Y[:, 0], label='HCE2')
ax.plot(states.t, states('LIGCC').Y[:, 0], label='LIGCC')
ax.plot(states.t, states('LIGOH').Y[:, 0], label='LIGOH')
ax.plot(states.t, states('LIG').Y[:, 0], label='LIG')
config(ax, xlabel='Time [s]', ylabel='Mass fraction [-]')

species = states.species_names
ys = [states(sp).Y[:, 0][-1] for sp in species]
ypos = np.arange(len(species))

fig, ax = plt.subplots(figsize=(6.4, 8), tight_layout=True)
ax.barh(ypos, ys, align='center')
ax.set_xlabel('Mass fraction [-]')
ax.set_ylim(min(ypos) - 1, max(ypos) + 1)
ax.set_yticks(ypos)
ax.set_yticklabels(species)
ax.invert_yaxis()
ax.set_axisbelow(True)
ax.set_frame_on(False)
ax.tick_params(color='0.8')
ax.xaxis.grid(True, color='0.8')

plt.show()

# fig, ax = plt.subplots(tight_layout=True)
# # ax.plot(states.t, states('CH2OHCHO').Y[:, 0], label='CH2OHCHO')
# # ax.plot(states.t, states('CHOCHO').Y[:, 0], label='CHOCHO')
# # ax.plot(states.t, states('CH3CHO').Y[:, 0], label='CH3CHO')
# # ax.plot(states.t, states('C6H6O3').Y[:, 0], label='C6H6O3')
# # ax.plot(states.t, states('C2H5CHO').Y[:, 0], label='C2H5CHO')
# config(ax, xlabel='Time [s]', ylabel='Mass fraction [-]')
