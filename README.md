## Kinetic Reactions for Biomass Pyrolysis
Various kinetic schemes are available to predict the devolatilization of biomass to gas, tar, and char products. The function files (`func_.py`) listed below are based on a particular kinetic scheme. Each function returns the product yield applicable to a particular scheme. See the comments in each file for more details. Different approaches of solving the system of kinetic reactions are also investigated in the `approach_.py` files.

*Requirements: Python 3, NumPy, and Matplotlib*

### funcChanBlasi.py
Kinetic functions based on Chan 1985 and Blasi 1993b kinetic reaction scheme for biomass pyrolysis. Primary and secondary reactions evaluated at some temperature. Separate function also provided just for primary reactions.

### References
* Chan, W.C.R., Kelbon, M. & Krieger, B.B., 1985. Modelling and experimental verification of physical and chemical processes during pyrolysis of a large biomass particle. Fuel, 64(11), pp.1505–1513.
* list
* list

### License
Files in this repository are available under the MIT license. See the LICENSE file for more info.
