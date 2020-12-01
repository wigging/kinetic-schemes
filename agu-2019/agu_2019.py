"""
Cantera batch reactor example where kinetic scheme based on Agu 2019 paper for
biomass fast pyrolysis.

Reference
---------
Cornelius E. Agu, Christoph Pfeifer, Marianne Eikeland, Lars-Andre Tokheim,
and Britt M.E. Moldestad. Detailed One-Dimensional Model for Steam-Biomass
Gasification in a Bubbling Fluidized Bed. Energy and Fuels, vol. 33, pp.
7385-7397, 2019.
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# Cantera batch reactor example with Blasi biomass pyrolysis kinetics
# ----------------------------------------------------------------------------

tk = 773.15                                 # temperature [K]
p = 101_325.0                               # pressure [Pa]
y = 'biomass:1 char:0 tar:0 volatiles:0'    # initial mass fractions [-]

gas = ct.Solution('agu.cti')
gas.TPY = tk, p, y
r = ct.IdealGasConstPressureReactor(gas, energy='off')

sim = ct.ReactorNet([r])
states = ct.SolutionArray(gas, extra=['t'])

# time vector to evaluate reaction rates [s]
time = np.linspace(0, 25, num=1000)

for tm in time:
    sim.advance(tm)
    states.append(r.thermo.state, t=tm)

# Print
# ----------------------------------------------------------------------------

print(f"""
Input parameters
----------------
T = {tk} K
P = {p} Pa""")

print("""
Final mass fractions
--------------------""")
for s in states.species_names:
    print(f'{s} = {states(s).Y[:, 0][-1]:.4f}')

# Plot
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states('biomass').Y[:, 0], label='biomass')
ax.plot(states.t, states('char').Y[:, 0], label='char')
ax.plot(states.t, states('tar').Y[:, 0], label='tar')
ax.plot(states.t, states('volatiles').Y[:, 0], label='volatiles')
ax.grid(color='0.9')
ax.legend(loc='best')
ax.set_frame_on(False)
ax.tick_params(color='0.9')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mass fraction [-]')

plt.show()
