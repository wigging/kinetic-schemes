# Kinetic Reaction Schemes for Biomass Pyrolysis

This repository contains various kinetic reaction schemes for biomass pyrolysis. Each file prepended with `kn` plots the products calculated from a particular scheme. Plots comparing all the schemes are provided by the `plot_primary` and `plot_secondary` files. Functions related to a particular scheme can be found in the `functions` folder. See the comments in each file and function for more details.

The scheme used in `kn_cpc_2016` is referred to as the "CPC scheme" which is the CPC's preferred kinetic scheme for woody biomass pyrolysis. It is a combination of the Blasi 1993, Chan 1985, and Liden 1988 kinetic schemes. The CPC scheme should be used as a starting point for computational models and as a basis of comparison between different models.

*Requirements: Python 3, NumPy, and Matplotlib*

## References

Kinetic schemes from the references listed below have been included in this repository. Files for a particular scheme are named according to author and year. For example, the scheme from the Di Blasi 1993 paper is plotted by the `kn_blasi_1993.py` file.

- Blasi, 1993. Combustion Science and Technology, 90, pp 315–340.

- Blasi, Branca, 2001. Ind. Eng. Chem. Res., 40, pp 5547-5556.

- Chan, Kelbon, Krieger, 1985. Fuel, 64, pp 1505–1513.

- Font, Marcilla, Verdu, Devesa, 1990. Ind. Egas. Chem. Res., 29, pp 1846-1855.

- Janse, Westerhout, Prins, 2000. Chem. Eng. Process., 39, pp 239-252.

- Koufopanos, 1991. The Canadian Journal of Chemical Engineering, 69, pp 907–915.

- Liden, Berruti, Scott, 1988. Chem. Eng. Comm., 65, pp 207-221.

- Miller, Bellan, 1997. Combust. Sci. and Tech., 126, pp 97-137.

- Papadikis, Gu, Bridgwater, 2010. Fuel Processing Technology, 91, pp 68–79.

- Ranzi, Corbetta, Pierucci, 2014. Chemical Engineering Science, 110, pp 2-12.

- Sadhukhan, Gupta, Saha, 2009. Bioresource Technology, 100, pp 3134-3139.

- Thurner, Mann, 1981. Ind. Eng. Chem. Process Des. Dev., 20, pp 482-488.

## License

Files in this repository are available under the MIT license. See the LICENSE file for more information.
