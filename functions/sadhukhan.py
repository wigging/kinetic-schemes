"""
Plot yields from primary and secondary reactions as determined by the
Sadhukhan 2009 kinetic scheme for biomass pyrolysis.

Reference:
Sadhukhan, Gupta, Saha, 2009. Bioresource Technology, 100, pp 3134-3139.
"""

import numpy as np

def sadhukhan(B, VG1, C1, VG2, C2, T, dt, s=1):
    """
    Primary and secondary kinetic reactions from Sadhukhan 2009 paper. Notice
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
    A1 = 168.4;     E1 = 51.965     # biomass -> (vol+gas)
    A2 = 13.2;      E2 = 45.960     # biomass -> char
    A3 = 5.7e6;     E3 = 92.4       # (vol+gas) + char -> (vol+gas)2 + char2
    R = 0.008314    # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R*T))    # biomass -> (volatiles + gases)1
    K2 = A2 * np.exp(-E2 / (R*T))    # biomass -> char1
    K3 = A3 * np.exp(-E3 / (R*T))    # (vol+gases)1 -> (vol+gases)2 + char2

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
