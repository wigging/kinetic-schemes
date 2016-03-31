"""
Function for Blasi/Chan/Liden kinetic scheme for biomass pyrolysis. This combined
scheme is referred to as the Cpc 2016 kinetic scheme. A similar scheme but
without water reaction was proposed in Papadikis 2010 paper. See comments in
function for more details.

References:
Blasi, 1993. Combustion Science and Technology, 90, pp 315–340.
Chan, Kelbon, Krieger, 1985. Fuel, 64, pp 1505–1513.
Liden, Berruti, Scott, 1988. Chem. Eng. Comm., 65, pp 207-221.
Papadikis, Gu, Bridgwater, 2010. Fuel Processing Technology, 91, pp 68–79.
"""

import numpy as np

def cpc(wood, gas, tar, char, water, vapor, T, dt, s=1):
    """
    Primary and secondary kinetic reactions for Cpc 2016 scheme based on
    Blasi 1993, Chan 1985, and Liden 1988 kinetics.

    Parameters
    ----------
    wood = wood concentration, kg/m^3
    gas = gas concentation, kg/m^3
    tar = tar concentation, kg/m^3
    char = char concentation, kg/m^3
    water = water concentration based on moisture content, kg/m^3
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
    nwater = new water concentration, kg/m^3
    nvapor = new water vapor concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.3e8;     E1 = 140    # wood -> gas   from Chan 1985
    A2 = 2.0e8;     E2 = 133    # wood -> tar   from Chan 1985
    A3 = 1.08e7;    E3 = 121    # wood -> char  from Chan 1985
    A4 = 4.28e6;    E4 = 107.5  # tar -> gas    from Liden 1988
    A5 = 1.0e6;     E5 = 108    # tar -> char   from Blasi 1993
    Aw = 5.13e6;    Ew = 87.9   # water -> water vapor  from Chan 1985
    R = 0.008314    # universal gas constant, kJ/mol*K

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))   # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))   # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))   # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T))   # tar -> gas
    K5 = A5 * np.exp(-E5 / (R * T))   # tar -> char
    Kw = Aw * np.exp(-Ew / (R * T))   # water -> vapor

    if s == 1:
        # primary reactions only
        rw = -(K1+K2+K3)*wood   # wood rate
        rg = K1*wood            # gas rate
        rt = K2*wood            # tar rate
        rc = K3*wood            # char rate
        rwt = -Kw*water         # moisture content rate
        rwv = Kw*water          # water vapor rate
        nwood = wood + rw*dt        # update wood concentration
        ngas = gas + rg*dt          # update gas concentration
        ntar = tar + rt*dt          # update tar concentration
        nchar = char + rc*dt        # update char concentration
        nwater = water + rwt*dt     # update water concentration
        nvapor = vapor + rwv*dt     # update water vapor concentation
    elif s == 2:
        # primary and secondary reactions
        rw = -(K1+K2+K3)*wood       # wood rate
        rg = K1*wood + K4*tar       # gas rate
        rt = K2*wood - (K4+K5)*tar  # tar rate
        rc = K3*wood + K5*tar       # char rate
        rwt = -Kw*water             # moisture content rate
        rwv = Kw*water              # water vapor rate
        nwood = wood + rw*dt            # update wood concentration
        ngas = gas + rg*dt              # update gas concentration
        ntar = tar + rt*dt              # update tar concentration
        nchar = char + rc*dt            # update char concentration
        nwater = water + rwt*dt         # update water concentration
        nvapor = vapor + rwv*dt         # update water vapor concentation

    # return new mass concentrations of products, kg/m^3
    return nwood, ngas, ntar, nchar, nwater, nvapor
