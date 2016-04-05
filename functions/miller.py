"""
Functions for Miller and Bellan 1997 kinetic scheme for biomass pyrolysis.
Biomass composition by mass of beech wood is used from Table 2 in paper. See
comments in each function for more details.

Reference:
Miller, Bellan, 1997. Combust. Sci. and Tech., 126, pp 97-137.
"""

import numpy as np

def millercell(cell, cella, gas, tar, char, T, dt, s=1):
    """
    Cellulose reactions from Miller and Bellan 1997 paper for biomass pyrolysis.

    Parameters
    ----------
    cell = cellulose concentration, kg/m^3
    cella = activie cellulose concentration, kg/m^3
    gas = gas concentration, kg/m^3
    tar = tar concentration, kg/m^3
    char = char concentration, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    ncell = new cellulose concentration, kg/m^3
    ncella = new active cellulose concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 2.8e19;    E1 = 242.4    # cell -> cella
    A2 = 3.28e14;   E2 = 196.5    # cella -> tar
    A3 = 1.3e10;    E3 = 150.5    # cella -> x*char + (1-x)*gas
    A4 = 4.28e6;    E4 = 108      # tar -> gas
    R = 0.008314    # universal gas constant, kJ/mol*K
    x = 0.35        # char formation mass ratio for cellulose

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # cell -> cella
    K2 = A2 * np.exp(-E2 / (R * T))  # cella -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # cella -> x*char + (1-x)*gas
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas

    if s == 1:
        # primary reactions only
        rcell = -K1*cell                    # cellulose rate
        rcella = K1*cell - (K2+K3)*cella    # active cellulose rate
        rgas = K3*cella                     # gas rate
        rtar = K2*cella                     # tar rate
        rchar = K3*cella                    # char rate
        ncell = cell + rcell*dt         # update cellulose concentration
        ncella = cella + rcella*dt      # update active cellulose concentration
        ngas = gas + rgas*(1-x)*dt      # update gas concentration
        ntar = tar + rtar*dt            # update tar concentration
        nchar = char + rchar*x*dt       # update char concentration
    elif s == 2:
        # primary and secondary reactions
        rcell = -K1*cell                    # cellulose rate
        rcella = K1*cell - (K2+K3)*cella    # active cellulose rate
        rgas = K3*cella + K4*tar            # gas rate
        rtar = K2*cella - K4*tar            # tar rate
        rchar = K3*cella                    # char rate
        ncell = cell + rcell*dt         # update cellulose concentration
        ncella = cella + rcella*dt      # update active cellulose concentration
        ngas = gas + rgas*(1-x)*dt      # update gas concentration
        ntar = tar + rtar*dt            # update tar concentration
        nchar = char + rchar*x*dt       # update char concentration

    # return new cellulose, active cellulose, gas, tar, char mass concentrations
    return ncell, ncella, ngas, ntar, nchar


def millerhemi(hemi, hemia, gas, tar, char, T, dt, s=1):
    """
    Hemicellulose reactions from Miller and Bellan 1997 paper for biomass
    pyrolysis.

    Parameters
    ----------
    hemi = hemicellulose concentration, kg/m^3
    hemia = active hemicellulose concentration, kg/m^3
    gas = gas concentration, kg/m^3
    tar = tar concentration, kg/m^3
    char = char concentration, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nhemi = new hemicellulose concentration, kg/m^3
    nhemia = new active hemicellulose concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 2.1e16;    E1 = 186.7    # hemi -> hemia
    A2 = 8.75e15;   E2 = 202.4    # hemia -> tar
    A3 = 2.6e11;    E3 = 145.7    # hemia -> x*char + (1-x)*gas
    A4 = 4.28e6;    E4 = 108      # tar -> gas
    R = 0.008314    # universal gas constant, kJ/mol*K
    x = 0.60        # char formation mass ratio for hemicellulose

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # hemi -> hemia
    K2 = A2 * np.exp(-E2 / (R * T))  # hemia -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # hemia -> x*char + (1-x)*gas
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas

    if s == 1:
        # primary reactions only
        rhemi = -K1*hemi                    # hemicellulose rate
        rhemia = K1*hemi - (K2+K3)*hemia    # active hemicellulose rate
        rgas = K3*hemia                     # gas rate
        rtar = K2*hemia                     # tar rate
        rchar = K3*hemia                    # char rate
        nhemi = hemi + rhemi*dt         # update hemicellulose concentration
        nhemia = hemia + rhemia*dt      # update active hemicellulose concentration
        ngas = gas + rgas*(1-x)*dt      # update gas concentration
        ntar = tar + rtar*dt            # update tar concentration
        nchar = char + rchar*x*dt       # update char concentration
    elif s == 2:
        # primary and secondary reactions
        rhemi = -K1*hemi                    # hemicellulose rate
        rhemia = K1*hemi - (K2+K3)*hemia    # active hemicellulose rate
        rgas = K3*hemia + K4*tar            # gas rate
        rtar = K2*hemia - K4*tar            # tar rate
        rchar = K3*hemia                    # char rate
        nhemi = hemi + rhemi*dt         # update hemicellulose concentration
        nhemia = hemia + rhemia*dt      # update active hemicellulose concentration
        ngas = gas + rgas*(1-x)*dt      # update gas concentration
        ntar = tar + rtar*dt            # update tar concentration
        nchar = char + rchar*x*dt       # update char concentration

    # return new hemicellulose, active hemicellulose, gas, tar, char mass concentrations
    return nhemi, nhemia, ngas, ntar, nchar


def millerlig(lig, liga, gas, tar, char, T, dt, s=1):
    """
    Lignin reactions from Miller and Bellan 1997 paper for biomass pyrolysis.

    Parameters
    ----------
    lig = lignin concentration, kg/m^3
    liga = active lignin concentration, kg/m^3
    gas = gas concentration, kg/m^3
    tar = tar concentration, kg/m^3
    char = char concentration, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nlig = new lignin concentration, kg/m^3
    nliga = new active lignin concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    nchar = new char concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 9.6e8;     E1 = 107.6      # lig -> liga
    A2 = 1.5e9;     E2 = 143.8      # liga -> tar
    A3 = 7.7e6;     E3 = 111.4      # liga -> x*char + (1-x)*gas
    A4 = 4.28e6;    E4 = 108        # tar -> gas
    R = 0.008314    # universal gas constant, kJ/mol*K
    x = 0.75        # char formation mass ratio for lignin

    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # lig -> liga
    K2 = A2 * np.exp(-E2 / (R * T))  # liga -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # liga -> x*char + (1-x)*gas
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas

    if s == 1:
        # primary reactions only
        rlig = -K1*lig                  # lignin rate
        rliga = K1*lig - (K2+K3)*liga   # active lignin rate
        rgas = K3*liga                  # gas rate
        rtar = K2*liga                  # tar rate
        rchar = K3*liga                 # char rate
        nlig = lig + rlig*dt        # update lignin concentration
        nliga = liga + rliga*dt     # update active lignin concentration
        ngas = gas + rgas*(1-x)*dt  # update gas concentration
        ntar = tar + rtar*dt        # update tar concentration
        nchar = char + rchar*x*dt   # update char concentration
    elif s == 2:
        # primary and secondary reactions
        rlig = -K1*lig                  # lignin rate
        rliga = K1*lig - (K2+K3)*liga   # active lignin rate
        rgas = K3*liga + K4*tar         # gas rate
        rtar = K2*liga - K4*tar         # tar rate
        rchar = K3*liga                 # char rate
        nlig = lig + rlig*dt        # update lignin concentration
        nliga = liga + rliga*dt     # update active lignin concentration
        ngas = gas + rgas*(1-x)*dt  # update gas concentration
        ntar = tar + rtar*dt        # update tar concentration
        nchar = char + rchar*x*dt   # update char concentration

    # return vector for each concentration, kg/m^3
    return nlig, nliga, ngas, ntar, nchar
