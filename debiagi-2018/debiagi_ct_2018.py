"""
Batch reactor example with Debiagi 2018 biomass pyrolysis kinetics.

This example uses Cantera to calculate conversion in a batch reactor
environment. All reactions and species are defined in the `cti` file where
`debiagi_sw.cti` is for softwood, `debiagi_hw.cti` is for hardwood, and
finally `debiagi_gr.cti` is for grass. Examples of initial mass fractions for
types of biomass such as softwood, hardwood, and grass are given below. Don't
forget to uncomment the relevant lines depending on the type of biomass.

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

# initial mass fractions for softwood [-] and specify cti file
# y0 = 'CELL:0.4002 GMSW:0.3112 LIGC:0.0287 LIGH:0.0196 LIGO:0.1379 TANN:0.0083 TGL:0.0941'
# gas = ct.Solution('debiagi_sw.cti')

# initial mass fractions for hardwood [-] and specify cti file
y0 = 'CELL:0.4109 XYHW:0.3150 LIGC:0.0202 LIGH:0.0133 LIGO:0.1443 TANN:0.0176 TGL:0.0786'
gas = ct.Solution('debiagi_hw.cti')

# initial mass fractions for grass [-] and specify cti file
# y0 = 'CELL:0.4232 XYGR:0.3465 LIGC:0.0059 LIGH:0.0167 LIGO:0.1086 TANN:0.0095 TGL:0.0845'
# gas = ct.Solution('debiagi_gr.cti')

# time vector to evaluate reaction rates [s]
time = np.linspace(0, 2.0, 100)

# Cantera batch reactor
# ----------------------------------------------------------------------------

gas.TPY = tk, p, y0
r = ct.IdealGasReactor(gas, energy='off')

sim = ct.ReactorNet([r])
states = ct.SolutionArray(gas, extra=['t'])

for tm in time:
    sim.advance(tm)
    states.append(r.thermo.state, t=tm)

# Print
# ----------------------------------------------------------------------------

print('--- Final mass fractions ---')
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
# ax.plot(states.t, states('GMSW').Y[:, 0], label='GMSW')
ax.plot(states.t, states('XYHW').Y[:, 0], label='XYHW')
# ax.plot(states.t, states('XYGR').Y[:, 0], label='XYGR')
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
