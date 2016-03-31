"""
Function for Koufopanos 1991 kinetic scheme for biomass pyrolysis. See comments
in function for more details.

Reference:
Koufopanos, 1991. The Canadian Journal of Chemical Engineering, 69, pp 907–915.
"""

import numpy as np

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
