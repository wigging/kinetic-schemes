"""
Function for the Liden 1988 kinetic scheme for biomass pyrolysis. See the
comments in the function for more details.

Reference:
Liden, Berruti, Scott, 1988. Chem. Eng. Comm., 65, pp 207-221.
"""

import numpy as np

def liden(wood, gas, tar, gaschar, T, dt, s=1):
    """
    Primary and secondary kinetic reactions from Liden 1988 paper. Parameters
    for total wood converstion (K = K1+K3) and secondary tar conversion (K2) are
    the only ones provided in paper. Can calculate K1 and K3 from phistar.

    Parameters
    ----------
    wood = wood concentration, kg/m^3
    gas = gas concentration, kg/m^3
    tar = tar concentation, kg/m^3
    gaschar = (gas+char) concentration, kg/m^3
    T = temperature, K
    dt = time step, s
    s = 1 primary reactions only, 2 primary and secondary reactions

    Returns
    -------
    nwood = new wood concentration, kg/m^3
    ngas = new gas concentration, kg/m^3
    ntar = new tar concentration, kg/m^3
    ngaschar = new (gas+char) concentration, kg/m^3
    """
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A = 1.0e13;         E = 183.3     # wood -> tar and (gas + char)
    A2 = 4.28e6;        E2 = 107.5    # tar -> gas
    R = 0.008314        # universal gas constant, kJ/mol*K
    phistar = 0.703     # maximum theoretical tar yield, (-)

    # reaction rate constant for each reaction, 1/s
    K = A * np.exp(-E / (R * T))        # wood -> tar and (gas + char)
    K1 = K * phistar                    # from phistar = K1/K
    K2 = A2 * np.exp(-E2 / (R * T))     # tar -> gas
    K3 = K - K1                         # from K = K1 + K3

    if s == 1:
        # primary reactions only
        rw = -K*wood        # wood rate
        rt = K1*wood        # tar rate
        rgc = K3*wood       # (gas+char) rate
        nwood = wood + rw*dt            # update wood concentration
        ngas = 0                        # no gas yield from primary reactions
        ntar = tar + rt*dt              # update tar concentration
        ngaschar = gaschar + rgc*dt     # update (gas+char) concentration
    elif s == 2:
        # primary and secondary reactions
        rw = -K*wood            # wood rate
        rg = K2*tar             # gas rate
        rt = K1*wood - K2*tar   # tar rate
        rgc = K3*wood           # (gas+char) rate
        nwood = wood + rw*dt            # update wood concentration
        ngas = gas + rg*dt              # update gas concentration
        ntar = tar + rt*dt              # update tar concentration
        ngaschar = gaschar + rgc*dt     # update (gas+char) concentration

    # return new wood, gas, tar, (gas+char) mass concentrations, kg/m^3
    return nwood, ngas, ntar, ngaschar
