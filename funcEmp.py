"""
Function to return empirical formula of biomass based on ultimate analysis.
"""

# Function 
# -----------------------------------------------------------------------------

def emp(c, h, o, n):
    """
    C = weight percent of carbon, %
    H = weight percent of hydrogen, %
    O = weight percent of oxygen, %
    N = weight percent of nitrogen, %
    note that the % should be entered as 48.1 not 0.481
    """
    
    # assume 100 grams of the compound, therefore wt. % = mass 
    mC = c
    mH = h
    mO = o
    mN = n
    
    # moles of C, H, O, N
    molC = mC/12.01
    molH = mH/1.01
    molO = mO/16.00
    molN = mN/14.01
    
    #return moles of C, H, O, N
    return molC, molH, molO, molN
    