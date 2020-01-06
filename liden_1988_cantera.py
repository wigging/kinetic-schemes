"""
Cantera batch reactor example where kinetic scheme based on Liden 1988 paper
for biomass fast pyrolysis.

Reference
---------
A.G. Liden, F. Berruti and D.S. Scott. A Kinetic Model for the Production of
Liquids from the Flash Pyrolysis of Biomass. Chemical Engineering
Communications, vol. 65, pp. 207-221, 1988.
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# Cantera batch reactor example with Liden 1988 biomass pyrolysis kinetics
# ----------------------------------------------------------------------------

tk = 773.15                             # temperature [K]
p = 101_325.0                           # pressure [Pa]
y = 'wood:1 gas:0 tar:0 gaschar:0'      # initial mass fractions [-]

gas = ct.Solution('liden.cti')
gas.TPY = tk, p, y
r = ct.IdealGasReactor(gas, energy='off')

sim = ct.ReactorNet([r])
states = ct.SolutionArray(gas, extra=['t'])

# time vector to evaluate reaction rates [s]
time = np.linspace(0, 25, num=1000)

for tm in time:
    sim.advance(tm)
    states.append(r.thermo.state, t=tm)

# assume a fixed carbon to estimate individual char and gas yields
phistar = 0.703             # maximum theoretical tar yield [-]
fc = 0.14                   # weight fraction of fixed carbon in wood
c3 = fc / (1 - phistar)     # char fraction for wood -> (gas+char)
g3 = 1 - c3                 # gas fraction for wood -> (gas+char)

# mass fractions for each component
wood = states('wood').Y[:, 0]
gas = states('gas').Y[:, 0] + g3 * states('gaschar').Y[:, 0]
tar = states('tar').Y[:, 0]
char = c3 * states('gaschar').Y[:, 0]

# total mass fractions should equal one
tot = wood + gas + tar + char

# Print
# ----------------------------------------------------------------------------

print(f"""
Input Parameters
----------------
T = {tk} K
P = {p} Pa
tmax = {25} s""")

print(f"""
Final Yields
------------""")

for s in states.species_names:
    print(f'{s} = {states(s).Y[:, 0][-1]:.4f}')

# Plot
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states('wood').Y[:, 0], label='wood')
ax.plot(states.t, states('gas').Y[:, 0], label='gas')
ax.plot(states.t, states('tar').Y[:, 0], label='tar')
ax.plot(states.t, states('gaschar').Y[:, 0], label='gaschar')
ax.grid(color='0.9')
ax.legend(loc='best')
ax.set_frame_on(False)
ax.tick_params(color='0.9')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mass fraction [-]')

fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, wood, label='wood')
ax.plot(states.t, gas, label='gas')
ax.plot(states.t, tar, label='tar')
ax.plot(states.t, char, label='char')
ax.plot(states.t, tot, 'k:', label='total')
ax.grid(color='0.9')
ax.legend(loc='best')
ax.set_frame_on(False)
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mass fraction [-]')
ax.tick_params(color='0.9')

plt.show()
