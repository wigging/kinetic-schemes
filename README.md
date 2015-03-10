## Kinetic Reactions for Biomass Pyrolysis
Various kinetic schemes are available to predict the devolatilization of biomass to gas, tar, and char products. The function files,`func_.py`, listed below are based on a particular kinetic scheme. Concentrations are on a volume basis represented as a density (kg/m<sup>3</sup>). Each function returns the product yield applicable to a particular scheme. See the comments in each file for more details.

*Requirements: Python 3, NumPy, and Matplotlib*

### approach.py

The `approach` files demonstrate different techniques in Python for solving a system of kinetic reactions for biomass pyrolysis. `approach_1.py` uses the analytical solution of the wood concentration as a function of time. `approach_2.py`  is an Euler method based on the initial wood concentration and reaction rates. `approach_3.py` uses the SciPy odeint to solve a system kinetic reactions as ODEs. `approach_4.py` uses the SciPy ode solver for a system of ODEs using techniques such as the Runge-Kutta method.

### funcBlasi.py

Contains functions based on the Blasi 1993 kinetic reaction scheme. The `blasi1` function returns the products for the primary reactions only. The `blasi2` function returns the primary and secondary reaction products.

### funcChan.py

Several functions based on the Chan 1985 kinetic reaction scheme. The `chan1` function returns only the products from the primary reactions. The `chan2` function returns the primary reaction products neglecting the moisture content to water vapor reaction. The `chan3` function returns the primary and secondary reaction products. The `chan4` function returns the primary and secondary reaction products without the moisture content reaction.

### funcFont.py

Two functions based on the Font 1990 kinetic scheme. The `font1` function represents the kinetic reactions developed from a fluidized bed experiment. The `font2` function represents the kinetic reactions developed from a pyroprobe experiment. The functions provide a comparison of how different experimental techniques affect the kinetic parameters.

### funcJanse.py

Functions based on the Janse 2000 kinetic reaction scheme. The `janse1` function returns the products for the primary reacions only. The `janse2` function returns the products for the primary and secondary reactions.

### funcPapadikis.py

Functions based on the Papadikis 2010 kinetic scheme. The `papadikis1` function returns the products for the primary reacions only. The `papadikis2` function returns the products for the primary and secondary reactions.

### funcThurner.py

A function based on the Thurner 1981 kinetic reaction scheme that returns the products from the primary reactions.

### References
* Di Blasi, C. (1993). Analysis of convection and secondary reaction effects within porous solid fuels undergoing pyrolysis. Combustion Science and Technology, 90(5-6), 315-340.
* Chan, W.C.R., Kelbon, M., & Krieger, B.B. (1985). Modelling and experimental verification of physical and chemical processes during pyrolysis of a large biomass particle. Fuel, 64(11), 1505-1513.
* Font, R., Marcilla, A., Verdu, E., & Devesa, J. (1990). Analysis of convection and secondary reaction effects within porous solid fuels undergoing pyrolysis. Ind. Eng. Chem. Res, 29, 1846-1855.
* Janse, A.M.C., Westerhout, R.W.J., & Prins, W. (2000). Modelling of flash pyrolysis of a single wood particle. Chemical Engineering and Processing: Process Intensification, 39(3), 239-252.
* Papadikis, K., Gu, S., & Bridgwater, A.V. (2010). Computational modelling of the impact of particle size to the heat transfer coefficient between biomass particles and a fluidised bed. Fuel Processing Technology, 91(1), 68-79.
* Thurner, F. & Mann, U., (1981). Kinetic Investigation of Wood Pyrolysis. Industrial and Engineering Chemistry Process Design and Development, 20, 482–488.

### License
Files in this repository are available under the MIT license. See the LICENSE file for more info.
