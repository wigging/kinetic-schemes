"""
Function for Thurner 1981 kinetic scheme for biomass pyrolysis. See comments in
the function for more details.

Reference:
Thurner, Mann, 1981. Ind. Eng. Chem. Process Des. Dev., 20, pp 482-488.
"""

import numpy as np

def thurner(wood, gas, tar, char, T, dt):
    """
    Primary kinetic reactions from Thurner 1981 paper for biomass pyrolysis.

    Parameters
    ----------
    wood = wood concentration, kg/m^3
    gas = gas concentation, kg/m^3
    tar = tar concentation, kg/m^3
    char = char concentation, kg/m^3
    T = temperature, K
    dt = time step, s

    Returns
    -------
    nwood = new wood concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.44e4;    E1 = 88.6     # wood -> gas
    A2 = 4.13e6;    E2 = 112.7    # wood -> tar
    A3 = 7.38e5;    E3 = 106.5    # wood -> char
    R = 0.008314    # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char

    # primary reactions only
    rw = -(K1+K2+K3) * wood     # wood rate
    rg = K1 * wood              # gas rate
    rt = K2 * wood              # tar rate
    rc = K3 * wood              # char rate
    nwood = wood + rw*dt    # update wood concentration
    ngas = gas + rg*dt      # update gas concentration
    ntar = tar + rt*dt      # update tar concentration
    nchar = char + rc*dt    # update char concentration

    # return new wood, char, gas, tar as mass concentrations, kg/m^3
    return nwood, ngas, ntar, nchar
