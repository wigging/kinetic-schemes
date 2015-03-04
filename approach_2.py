"""
Euler approach to solve system of kinetic reactions for biomass pyrolysis.
Solution based on initial wood density and reaction rates, for example
pw = pwi + rww*dt.

Requirements:
Pyton 3, Numpy, Matplotlib

References:
1) Papadikis, Gu, Bridgwater, 2010. Fuel Processing Technology, 91(1), pp.68–79.
2) Chan, Kelbon, Krieger, 1985. Fuel, 64(11), pp.1505–1513.
3) Liden, Berruti, Scott, 1988. Chemical Engineering Communications, 65, pp.207–221.
4) Blasi, 1993. Combustion Science and Technology, 90, pp.315–340.
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py

py.close('all')

# Euler function to solve for wood, gas, tar, char concentrations
#------------------------------------------------------------------------------

def euler(T, pw, pg, pt, pc, dt, i):
    """
    Function for Euler approach to solve a system of kinetic reactions as ODEs.
    INPUTS: 
    T = temperature, K                      pw = wood concentration, kg/m^3
    pg = gas concentration, kg/m^3          pt = tar concentration, kg/m^3
    pc = char concentration, kg/m^3         pm = moisture concentration, kg/m^3
    pv = water vapor concentration, kg/m^3  dt = delta time, s
    i = index    
    OUTPUTS:
    pww = wood  pgg = gas   ptt = tar   pcc = char  pmm = moisture
    pvv = water vapor
    """
    
    R = 0.008314 # universal gas constant, kJ/mol*K
    
    # Kinetic parameters from Chan 1985 and Blasi 1993b
    # A = pre-factor (1/s) and E = activation energy (kJ/mol)
    A1 = 1.3e8;     E1 = 140    # wood -> gas
    A2 = 2e8;       E2 = 133    # wood -> tar
    A3 = 1.08e7;    E3 = 121    # wood -> char
    A4 = 4.28e6;    E4 = 108    # tar -> gas
    A5 = 1e6;       E5 = 108    # tar -> char
    
    # reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp(-E1 / (R * T))  # wood -> gas
    K2 = A2 * np.exp(-E2 / (R * T))  # wood -> tar
    K3 = A3 * np.exp(-E3 / (R * T))  # wood -> char
    K4 = A4 * np.exp(-E4 / (R * T))  # tar -> gas
    K5 = A5 * np.exp(-E5 / (R * T))  # tar -> char
    
    # reaction rates where r = dp/dt, rho/s or kg/(m^3 s)
    rw = -(K1+K2+K3)*pw[i-1]                    # wood rw = dpw/dt
    rg = K1*pw[i-1] +  K4*pt[i-1]               # gas rg = dpg/dt
    rt = K2*pw[i-1] - K4*pt[i-1] - K5*pt[i-1]   # tar rt = dpt/dt
    rc = K3*pw[i-1] + K5*pt[i-1]                # char rc = dpc/dt
    
    # update concentrations as a density, kg/m^3
    pww = pw[i-1] + rw*dt   # wood pw
    pgg = pg[i-1] + rg*dt   # gas pg
    ptt = pt[i-1] + rt*dt   # tar pt
    pcc = pc[i-1] + rc*dt   # char pc
    
    # return new wood, gas, tar, char concentrations as a density, kg/m^3
    return pww, pgg, ptt, pcc


# Parameters from Papadikis 2010a
#------------------------------------------------------------------------------

rhow = 700  # density of wood, kg/m^3
Tinf = 773  # ambient temp, K

# Initial Calculations
#------------------------------------------------------------------------------

dt = 0.01   # time step, delta t
tmax = 25   # max time, s
t = np.linspace(0, tmax, num=tmax/dt)   # time vector
p = len(t)  # total number of time steps

# Calculate Kinetic Reactions and Concentrations
#------------------------------------------------------------------------------

# vectors for wood, gas, tar, char concentrations as a density, kg/m^3
pw = np.zeros(len(t))   # wood 
pg = np.zeros(len(t))   # gas
pt = np.zeros(len(t))   # tar
pc = np.zeros(len(t))   # char

pw[:] = rhow # initial wood density

# solve kinetics for primary and secondary reactions
for i in range(1, p):
    pw[i], pg[i], pt[i], pc[i] = euler(Tinf, pw, pg, pt, pc, dt, i)
 
# Plot Results
#------------------------------------------------------------------------------

py.rcParams['xtick.major.pad'] = 8
py.rcParams['ytick.major.pad'] = 8
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

py.figure(2)
py.plot(t, pw, label='wood')
py.plot(t, pg, label='gas')
py.plot(t, pt, label='tar')
py.plot(t, pc, label='char')
py.legend(loc='best', numpoints=1)
py.title('Euler Approach, reactions at T = %.f K' % Tinf)
py.xlabel('Time (s)')
py.ylabel('Concentration ($kg/m^3$)')
py.show()
