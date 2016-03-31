"""
Functions for Blasi 1993 and Blasi Branca 2001 kinetic reaction schemes for
biomass pyrolysis. See comments in each function for more details.

References:
Blasi, 1993. Combustion Science and Technology, 90, pp 315–340.
Blasi, Branca, 2001. Ind. Eng. Chem. Res., 40, pp 5547-5556.
"""

import numpy as np

def blasi(wood, gas, tar, char, T, dt, s=1):
    """
    Primary and secondary kinetic reactions from Table 1 in Blasi 1993 paper.
    Note that primary reaction parameters in table are not cited correctly from
    the Thurner and Mann 1981 paper, this function uses the correct parameters.

    Parameters
    ----------
    wood = wood concentration, kg/m^3
    gas = gas concentation, kg/m^3
    tar = tar concentation, kg/m^3
    char = char concentation, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nwood = new wood concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.4345e4;      E1 = 88.6     # wood -> gas
    A2 = 4.125e6;       E2 = 112.7    # wood -> tar
    A3 = 7.3766e5;      E3 = 106.5    # wood -> char
    A4 = 4.28e6;        E4 = 108      # tar -> gas
    A5 = 1.0e6;         E5 = 108      # tar -> char
    R = 0.008314        # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas
    K5 = A5 * np.exp(-E5 / (R * T))  # tar -> char

    if s == 1:
        # primary reactions only
        rw = -(K1+K2+K3)*wood   # wood rate
        rg = K1*wood            # gas rate
        rt = K2*wood            # tar rate
        rc = K3*wood            # char rate
        nwood = wood + rw*dt    # update wood concentration
        ngas = gas + rg*dt      # update gas concentration
        ntar = tar + rt*dt      # update tar concentration
        nchar = char + rc*dt    # update char concentration
    elif s == 2:
        # primary and secondary reactions
        rw = -(K1+K2+K3)*wood       # wood rate
        rg = K1*wood + K4*tar       # gas rate
        rt = K2*wood - (K4+K5)*tar  # tar rate
        rc = K3*wood + K5*tar       # char rate
        nwood = wood + rw*dt    # update wood concentration
        ngas = gas + rg*dt      # update gas concentration
        ntar = tar + rt*dt      # update tar concentration
        nchar = char + rc*dt    # update char concentration

    # return new wood, gas, tar, char mass concentrations, kg/m^3
    return nwood, ngas, ntar, nchar


def blasibranca(pw, pg, pt, pc, T, dt):
    """
    Primary kinetic reactions from Table 1 in Blasi and Branca 2001 paper.

    Parameters
    ----------
    pw = wood concentration, kg/m^3
    pg = gas concentation, kg/m^3
    pt = tar concentation, kg/m^3
    pc = char concentation, kg/m^3
    T = temperature, K
    dt = time step, s

    Returns
    -------
    nw = new wood concentration, kg/m^3
    ng = new gas concentration, kg/m^3
    nt = new tar concentration, kg/m^3
    nc = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 4.38e9;    E1 = 152.7    # wood -> gas
    A2 = 1.08e10;   E2 = 148      # wood -> tar
    A3 = 3.27e6;    E3 = 111.7    # wood -> char
    R = 0.008314 # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char

    # primary reactions only
    rw = -(K1+K2+K3)*pw # wood rate
    rg = K1*pw          # gas rate
    rt = K2*pw          # tar rate
    rc = K3*pw          # char rate
    nw = pw + rw*dt # update wood concentration
    ng = pg + rg*dt # update gas concentration
    nt = pt + rt*dt # update tar concentration
    nc = pc + rc*dt # update char concentration

    # return new wood, gas, tar, char as mass concentrations, kg/m^3
    return nw, ng, nt, nc
