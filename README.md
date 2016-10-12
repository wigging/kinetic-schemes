# Kinetic Reaction Schemes for Biomass Pyrolysis

This repository contains various kinetic reaction schemes for biomass
pyrolysis. Each file plots the products for that particular scheme. For
example, `blasi_2001.py` uses the kinetic parameters from the Di Blasi 2001
paper to plot the conversion of wood to gas, tar, and char products. Plots
comparing all the kinetic schemes are provided by the `plot_primary` and
`plot_secondary` files which use the functions contained in the `functions`
folder. See the comments in each file for more details.

The scheme used in `cpc_2016` is referred to as the "CPC scheme". It is a
combination of the Blasi 1993, Chan 1985, and Liden 1988 kinetic schemes. The
CPC scheme should be used as a starting point for computational models and as a
basis of comparison between different types of models.

Contributions to this repository can be made via a Pull Request or suggested as
an Issue. Our contact information is also available on the
[cpcbiomass.org](http://cpcbiomass.org) website.

Python 3, NumPy, and Matplotlib are needed to run the examples in this repo.

## References

Kinetic schemes from the references listed below have been included in this
repository. Files for a particular scheme are named according to author and
year of publication. For example, the scheme from the Di Blasi 1993 paper is plotted by the
`blasi_1993.py` file.

- Adjaye, 1993. Thesis Appendix F, pp 311-319.

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

Files in this repository are available under the MIT license. See the LICENSE
file for more information.

