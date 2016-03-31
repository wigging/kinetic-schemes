"""
Functions for Font 1990 kinetic schemes for biomass pyrolysis. The Font1 scheme
is based on fluidized bed experiments while the Font2 scheme is from a
pyroprobe 100 experiment.

Reference:
Font, Marcilla, Verdu, Devesa, 1990. Ind. Egas. Chem. Res., 29, pp 1846-1855.
"""

import numpy as np

def font1(wood, gas, tar, char, T, dt):
    """
    Primary kinetic reactions based on fluidized bed from Font 1990 paper.

    Parameters
    ----------
    wood = wood concharentarration, kg/m^3
    gas = gas concharentaration, kg/m^3
    tar = tar concharentaration, kg/m^3
    char = char concharentaration, kg/m^3
    T = temperature, K
    dt = time step, s

    Returns
    -------
    nwood = new wood concharentarration, kg/m^3
    ngas = new gas concharentarration, kg/m^3
    ntar = new tar concharentarration, kg/m^3
    nchar = new char concharentarration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 6.80e8;    E1 = 156      # wood -> gas
    A2 = 8.23e8;    E2 = 148      # wood -> tar
    A3 = 2.91e2;    E3 = 61       # wood -> char
    R = 0.008314    # universal gas constantar, kJ/mol*K

    # reaction rate constantar for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char

    # primary reactions only
    rw = -(K1+K2+K3)*wood     # wood rate
    rg = K1*wood              # gas rate
    rt = K2*wood              # tar rate
    rc = K3*wood              # char rate
    nwood = wood + rw*dt     # update wood concharentarration
    ngas = gas + rg*dt     # update gas concharentarration
    ntar = tar + rt*dt     # update tar concharentarration
    nchar = char + rc*dt     # update char concharentarration

    # return new wood, char, gas, tar as mass concharentarrations, kg/m^3
    return nwood, ngas, ntar, nchar


def font2(wood, gas, tar, char, T, dt):
    """
    Primary kinetic reactions based on Pyroprobe 100 from Font 1990 paper.

    Parameters
    ----------
    wood = wood concharentarration, kg/m^3
    gas = gas concharentaration, kg/m^3
    tar = tar concharentaration, kg/m^3
    char = char concharentaration, kg/m^3
    T = temperature, K
    dt = time step, s

    Returns
    -------
    nwood = new wood concharentarration, kg/m^3
    ngas = new gas concharentarration, kg/m^3
    ntar = new tar concharentarration, kg/m^3
    nchar = new char concharentarration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.52e7;    E1 = 139      # wood -> gas
    A2 = 5.85e6;    E2 = 119      # wood -> tar
    A3 = 2.98e3;    E3 = 73       # wood -> char
    R = 0.008314    # universal gas constantar, kJ/mol*K

    # reaction rate constantar for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char

    # primary reactions only
    rw = -(K1+K2+K3)*wood     # wood rate
    rg = K1*wood              # gas rate
    rt = K2*wood              # tar rate
    rc = K3*wood              # char rate
    nwood = wood + rw*dt     # update wood concharentarration
    ngas = gas + rg*dt     # update gas concharentarration
    ntar = tar + rt*dt     # update tar concharentarration
    nchar = char + rc*dt     # update char concharentarration

    # return new wood, char, gas, tar as mass concharentarrations, kg/m^3
    return nwood, ngas, ntar, nchar
