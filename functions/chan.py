"""
Function for Chan 1985 kinetic reaction scheme for biomass pyrolysis. See
comments in the function for more details.

Reference:
Chan, Kelbon, Krieger, 1985. Fuel, 64, pp 1505–1513.
"""

import numpy as np

def chan(wood, gas, tar, char, water, vapor, T, dt, s=1):
    """
    Primary and secondary reactions with moisture reaction from Table 2 in the
    Chan 1985 paper.

    Parameters
    ----------
    wood = wood concentration, kg/m^3
    gas = gas concentation, kg/m^3
    tar = tar concentation, kg/m^3
    char = char concentation, kg/m^3
    water = moisture content concentration, kg/m^3
    vapor = water vapor concentration, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nwood = new wood concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    nwater = new moisture concentration, kg/m^3
    nvapor = new water vapor concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.3e8;     E1 = 140    # wood -> gas1
    A2 = 2.0e8;     E2 = 133    # wood -> tar1
    A3 = 1.08e7;    E3 = 121    # wood -> char
    A4 = 5.13e6;    E4 = 87.9   # moisture -> water vapor
    A5 = 1.48e6;    E5 = 144    # tar1 -> gas2 + tar2
    R = 0.008314    # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas1
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar1
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T))  # moisture -> water vapor
    K5 = A5 * np.exp(-E5 / (R * T))  # tar1 -> gas2 + tar2

    if s == 1:
        # primary reactions only
        rw = -(K1+K2+K3)*wood     # wood rate
        rg = K1*wood              # gas rate
        rt = K2*wood              # tar rate
        rc = K3*wood              # char rate
        rm = -K4*water             # moisture content rate
        rv = K4*water              # water vapor rate
        nwood = wood + rw*dt     # update wood concentration
        ngas = gas + rg*dt     # update gas concentration
        ntar = tar + rt*dt     # update tar concentration
        nchar = char + rc*dt     # update char concentration
        nwater = water + rm*dt     # update moisture content concentration
        nvapor = vapor + rv*dt     # update water vapor concentration
    elif s == 2:
        # primary and secondary reactions
        rw = -(K1+K2+K3)*wood     # wood rate
        rg = K1*wood + K5*tar      # gas rate
        rt = K2*wood - K5*tar      # tar rate
        rc = K3*wood              # char rate
        rm = -K4*water             # moisture content rate
        rv = K4*water              # water vapor rate
        nwood = wood + rw*dt     # update wood concentration
        ngas = gas + rg*dt     # update gas concentration
        ntar = tar + rt*dt     # update tar concentration
        nchar = char + rc*dt     # update char concentration
        nwater = water + rm*dt     # update moisture content concentration
        nvapor = vapor + rv*dt     # update water vapor concentration

    # return wood, gas, tar, char, moisture, water vapor as mass concentrations
    return nwood, ngas, ntar, nchar, nwater, nvapor
