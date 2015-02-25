"""
Kinetic Reaction Scheme Functions for Fast Pyrolysis of Biomass.

Each function is for a particular kinetic scheme.
Reference for each scheme is provided as main author and publication year.
"""

# modules
# -----------------------------------------------------------------------------
import numpy as np

# Sadhukhan2009 
# volatiles+gases, char, primary and secondary reactions
# -----------------------------------------------------------------------------

def kn1(T, pw, pc, pg, dt, i, H):
    
    R = 0.008314 # universal gas constant, kJ/mol*K
    
    # A as pre-factor (1/s) and E as activation energy (kJ/mol)
    A1 = 168.4; E1 = 51.965     # biomass -> volatiles + gases
    A2 = 13.2;  E2 = 45.960     # biomass -> char
    A3 = 5.7e6; E3 = 92.4       # (vol+gases)1 -> (vol+gases)2
    
    # evaluate reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T[i]))  # biomass -> volatiles + gases
    K2 = A2 * np.exp(-E2 / (R * T[i]))  # biomass -> char
    K3 = A3 * np.exp(-E3 / (R * T[i]))  # (vol+gases)1 -> (vol+gases)2
    
    # determine reaction rate for each reaction, rho/s
    rw = -(K1+K2) * pw[i-1]                     # wood rate
    rg1 = K1 * pw[i-1] - K3*pg[i-1] * pc[i-1]   # gas 1 rate
    rc1 = K2 * pw[i-1] - K3*pg[i-1] * pc[i-1]   # char 1 rate
    rg2 = K3 * pg[i-1] * pc[i-1]                # gas 2 rate
    rc2 = K3 * pg[i-1] * pc[i-1]                # char 2 rate
    
    # update wood, char, gas concentration as a density, kg/m^3
    pww = pw[i-1] + rw * dt             # wood
    pcc = pc[i-1] + (rc1 + rc2) * dt    # char
    pgg = pg[i-1] + (rg1 + rg2) * dt    # gas
    
    # calculate heat of generation term
    rp = -K1*pww    # rate of pyrolysis
    g = H*rp        # heat generation
    
    # return the wood, char, gas concentration and the heat of generation
    return pww, pcc, pgg, g


# Chan1985, Blasi1993b
# primary and secondary reactions
# -----------------------------------------------------------------------------

def kn2(T, pw, pc, pg, pt, dt, i, H):
    
    R = 0.008314 # universal gas constant, kJ/mol*K
    
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.3e8;  E1 = 140    # wood -> gas
    A2 = 2e8;    E2 = 133    # wood -> tar
    A3 = 1.08e7; E3 = 121    # wood -> char
    A4 = 4.28e6; E4 = 108    # tar -> gas
    A5 = 1e6;    E5 = 108    # tar -> char
    
    # evaluate reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T[i]))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T[i]))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T[i]))  # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T[i]))  # tar -> gas
    K5 = A5 * np.exp(-E5 / (R * T[i]))  # tar -> char
    
    # determine reaction rate for each reaction, rho/s
    rww = -(K1+K2+K3) * pw[i-1]     # wood rate
    rwg = K1 * pw[i-1]              # wood -> gas rate
    rwt = K2 * pw[i-1]              # wood -> tar rate
    rwc = K3 * pw[i-1]              # wood -> char rate
    rtg = K4 * pt[i-1]              # tar -> gas rate
    rtc = K5 * pt[i-1]              # tar -> char rate
    
    # update wood, char, gas concentration as a density, kg/m^3
    pww = pw[i-1] + rww*dt                  # wood
    pgg = pg[i-1] + (rwg + rtg)*dt          # gas
    ptt = pt[i-1] + (rwt - rtg - rtc)*dt    # tar
    pcc = pc[i-1] + (rwc + rtc)*dt          # char
    
    # calculate heat of generation term
    g = H*rww        # heat generation, W/m^3
    
    # return the wood, char, gas, tar concentration and the heat generation
    return pww, pcc, pgg, ptt, g


# Chan1985
# moisture content, heat of vaporization, no secondary reactions
# -----------------------------------------------------------------------------

def kn3(T, pw, pc, pg, pt, pwa, pva, dt, i, H):
    
    R = 0.008314 # universal gas constant, kJ/mol*K
    
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.3e8;     E1 = 140    # wood -> gas
    A2 = 2e8;       E2 = 133    # wood -> tar
    A3 = 1.08e7;    E3 = 121    # wood -> char
    Aw = 5.13e6;    Ew = 87.9   # water -> vapor
    
    # evaluate reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T[i]))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T[i]))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T[i]))  # wood -> char
    Kw = Aw * np.exp(-Ew / (R * T[i]))  # water -> vapor
    
    # determine reaction rate for each reaction, rho/s
    rww = -(K1+K2+K3) * pw[i-1]     # rate of wood pyrolysis
    rwg = K1 * pw[i-1]              # rate of wood -> gas
    rwt = K2 * pw[i-1]              # rate of wood -> tar
    rwc = K3 * pw[i-1]              # rate of wood -> char
    rwa = -Kw * pwa[i-1]            # rate of water vaporization
    rva = Kw * pwa[i-1]             # rate of water -> vapor
    
    # update concentrations as a density, kg/m^3
    pww = pw[i-1] + rww*dt      # wood
    pgg = pg[i-1] + rwg*dt      # gas
    ptt = pt[i-1] + rwt*dt      # tar
    pcc = pc[i-1] + rwc*dt      # char
    pwwa = pwa[i-1] + rwa*dt    # water
    pvva = pva[i-1] + rva*dt    # vapor
    
    # calculate heat of generation term
    Hv = 2260000        # heat of vaporization, J/kg
    g = H*rww + Hv*rwa  # heat generation, W/m^3    
    
    # return wood, char, gas, tar, water, vapor concentration & heat generation
    return pww, pcc, pgg, ptt, pwwa, pvva, g
    
    
# Chan1985, Blasi1993b
# moisture content, heat of vaporization, primary and secondary reactions
# -----------------------------------------------------------------------------

def kn4(T, pw, pc, pg, pt, pwa, pva, dt, i, H):
    
    R = 0.008314 # universal gas constant, kJ/mol*K
    
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.3e8;     E1 = 140    # wood -> gas
    A2 = 2e8;       E2 = 133    # wood -> tar
    A3 = 1.08e7;    E3 = 121    # wood -> char
    A4 = 4.28e6;    E4 = 108    # tar -> gas
    A5 = 1e6;       E5 = 108    # tar -> char
    Aw = 5.13e6;    Ew = 87.9   # water -> vapor
    
    # evaluate reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T[i]))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T[i]))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T[i]))  # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T[i]))  # tar -> gas
    K5 = A5 * np.exp(-E5 / (R * T[i]))  # tar -> char
    Kw = Aw * np.exp(-Ew / (R * T[i]))  # water -> vapor
    
    # determine reaction rate for each reaction, rho/s
    rww = -(K1+K2+K3) * pw[i-1]     # wood rate
    rwg = K1 * pw[i-1]              # wood -> gas rate
    rwt = K2 * pw[i-1]              # wood -> tar rate
    rwc = K3 * pw[i-1]              # wood -> char rate
    rtg = K4 * pt[i-1]              # tar -> gas rate
    rtc = K5 * pt[i-1]              # tar -> char rate
    rwa = -Kw * pwa[i-1]            # rate of water vaporization
    rva = Kw * pwa[i-1]             # rate of water -> vapor
    
    # update wood, char, gas concentration as a density, kg/m^3
    pww = pw[i-1] + rww*dt                  # wood
    pgg = pg[i-1] + (rwg + rtg)*dt          # gas
    ptt = pt[i-1] + (rwt - rtg - rtc)*dt    # tar
    pcc = pc[i-1] + (rwc + rtc)*dt          # char
    pwwa = pwa[i-1] + rwa*dt                # water
    pvva = pva[i-1] + rva*dt                # vapor
    
    # calculate heat of generation term
    Hv = 2260000        # heat of vaporization, J/kg
    g = H*rww + Hv*rwa  # heat generation, W/m^3  
    
    # return the wood, char, gas, tar concentration and the heat generation
    return pww, pcc, pgg, ptt, pwwa, pvva, g
    
    