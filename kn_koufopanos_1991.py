"""
Plot yields from primary and secondary reactions as determined by the
Koufopanos 1991 kinetic scheme for biomass pyrolysis. Note that this scheme
focuses on wood conversion and char yield. Product of volatiles and gas is
lumped together as (V+G) so individual tar and gas component is not provided.

Reference:
Koufopanos, 1991. The Canadian Journal of Chemical Engineering, 69, pp 907–915.
"""

import numpy as np
import matplotlib.pyplot as py

# Parameters
# ------------------------------------------------------------------------------

T = 773  # temperature for rate constants, K

dt = 0.01                               # time step, delta t
tmax = 25                               # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
nt = len(t)                             # total number of time steps

# Function for Koufopanos 1991 Kinetic Scheme
# ------------------------------------------------------------------------------

def koufopanos(B, VG1, C1, VG2, C2, T, dt, s=1):
    """
    Primary and secondary kinetic reactions from Koufopanos 1991 paper. Notice
    that volatiles and gases are grouped together as (Volatiles + Gases) which
    are labeled here as VG.

    Parameters
    ----------
    B = biomass concentration
    VG1 = (volatiles + gas)1 concentration
    C1 = char1 concentration
    VG2 = (volatiles + gas)2 concentration
    C2 = char2 concentration
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nB = new biomass concentration
    nVG1 = new (volatiles + char)1 concentration
    nC1 = new char1 concentration
    nVG2 = new (volatiles + gas)2 concentration
    nC2 = new char2 concentration
    """
    # A as pre-factor (1/s) and E as activation energy (kJ/mol)
    A1 = 9.973e-5;  G1 = 17254.4;   L1 = -9061227   # biomass -> (volatiles + gases)1
    A2 = 1.068e-3;  G2 = 10224.4;   L2 = -6123081   # biomass -> char1
    A3 = 5.7e5;     E3 = 81         # (vol+gases)1 -> (vol+gases)2 + char2
    R = 0.008314    # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp((G1 / T) + (L1 / T**2))    # biomass -> (volatiles + gases)1
    K2 = A2 * np.exp((G2 / T) + (L2 / T**2))    # biomass -> char1
    K3 = A3 * np.exp(-E3 / (R * T))             # (vol+gases)1 -> (vol+gases)2 + char2

    if s == 1:
        # primary reactions only
        rB = -(K1+K2)*B     # biomass rate
        rVG1 = K1*B         # (volatiles + gases)1 rate
        rC1 = K2*B          # char1 rate
        rVG2 = 0            # (volatiles + gases)2 rate
        rC2 = 0             # char2 rate
        nB = B + rB*dt          # update biomass concentration
        nVG1 = VG1 + rVG1*dt    # update (volatiles + gases)1 concentration
        nC1 = C1 + rC1*dt       # update char1 concentration
        nVG2 = VG2 + rVG2*dt    # update (volatiles + gases)2 concentration
        nC2 = C2 + rC2*dt       # update char2 concentration
    elif s == 2:
        # primary and secondary reactions
        rB = -(K1+K2)*B     # biomass rate
        rVG1 = K1*B         # volatiles + gases)1 rate
        rC1 = K2*B - K3*C1  # char1 rate
        rVG2 = K3*C1        # (volatiles + gases)2 rate
        rC2 = K3*C1         # char2 rate
        nB = B + rB*dt          # update biomass concentration
        nVG1 = VG1 + rVG1*dt    # update (volatiles + gases)1 concentration
        nC1 = C1 + rC1*dt       # update char1 concentration
        nVG2 = VG2 + rVG2*dt    # update (volatiles + gases)2 concentration
        nC2 = C2 + rC2*dt       # update char2 concentration

    return nB, nVG1, nC1, nVG2, nC2

# Product from Kinetic Scheme
# ------------------------------------------------------------------------------

# Assume initial concentration of B(0) = 1 and everything else initially at zero
# such as VG(0) = C(0) = 0 where VG is (Volatiles + Gases) and C is Char.

# concentrations reported on a mass basis in kg/m^3
# pw = wood concentration, pg = gas concentration, pt = tar concentration,
# pc = char concentration

# store concentrations from primary reactions at each time step
B = np.ones(nt)   # biomass concentration vector
VG1 = np.zeros(nt)     # (volatiles + gases)1 concentration vector
C1 = np.zeros(nt)      # char1 concentration vector
VG2 = np.zeros(nt)     # (volatiles + gases)2 concentration vector
C2 = np.zeros(nt)      # char2 concentration vector

# store concentrations from primary and secondary reactions at each time step
B_2 = np.ones(nt)   # biomass concentration vector
VG1_2 = np.zeros(nt)     # (volatiles + gases)1 concentration vector
C1_2 = np.zeros(nt)      # char1 concentration vector
VG2_2 = np.zeros(nt)     # (volatiles + gases)2 concentration vector
C2_2 = np.zeros(nt)      # char2 concentration vector

# products from primary reactions only
for i in range(1, nt):
    B[i], VG1[i], C1[i], VG2[i], C2[i] = koufopanos(B[i-1], VG1[i-1], C1[i-1], VG2[i-1], C2[i-1], T, dt)

# products from primary and secondary reactions
for i in range(1, nt):
    B_2[i], VG1_2[i], C1_2[i], VG2_2[i], C2_2[i] = koufopanos(B_2[i-1], VG1_2[i-1], C1_2[i-1], VG2_2[i-1], C2_2[i-1], T, dt, s=2)

# totals from primary reactions only
pvg = VG1 + VG2
pc = C1 + C2

# totals from primary and secondary reactions, assume VG1 -> (VG + C)2 where
# components in the group (VG + C)2 = 1/2*VG2 + 1/2*C2
pvg_2 = VG1_2 + 0.5*VG2_2
pc_2 = C1_2 + 0.5*C2_2

# mass balance to check results
mt = B + pvg + pc
mt2 = B_2 + pvg_2 + pc_2

# Plot Results
# ------------------------------------------------------------------------------

py.ion()
py.close('all')

py.figure(1)
py.plot(t, B, lw=2, label='B')
py.plot(t, pvg, lw=2, label='(V+G)$_1$')
py.plot(t, pc, lw=2, label='Char$_1$')
py.title('Koufopanos 1991 primary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration (normalized mass basis)')
py.legend(loc='best', numpoints=1)
py.grid()

py.figure(2)
py.plot(t, B_2, lw=2, label='B')
py.plot(t, pvg_2, lw=2, label='(V+G)')
py.plot(t, pc_2, lw=2, label='Char')
py.title('Koufopanos 1991 primary and secondary reactions at T = {} K'.format(T))
py.xlabel('Time (s)')
py.ylabel('Concentration (normalized mass basis)')
py.legend(loc='best', numpoints=1)
py.grid()
