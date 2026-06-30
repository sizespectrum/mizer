# Calculate the number of individuals within a size range

Calculates the number of individuals within user-defined size limits.
The default option is to use the whole size range. You can specify
minimum and maximum weight or lengths for the species. Lengths take
precedence over weights (i.e. if both min_l and min_w are supplied, only
min_l will be used)

## Usage

``` r
getN(object, ...)
```

## Arguments

- object:

  An object of class `MizerParams` or `MizerSim`.

- ...:

  Arguments passed on to
  [`get_size_range_array`](https://sizespectrum.org/mizer/reference/get_size_range_array.md)

  `min_w`

  :   Smallest weight in size range. Defaults to smallest weight in the
      model.

  `max_w`

  :   Largest weight in size range. Defaults to largest weight in the
      model.

  `min_l`

  :   Smallest length in size range. If supplied, this takes precedence
      over `min_w`.

  `max_l`

  :   Largest length in size range. If supplied, this takes precedence
      over `max_w`.

## Value

If called with a MizerParams object, a named vector with the numbers for
each species in the model. If called with a MizerSim object, a
`ArrayTimeBySpecies` object (time x species) containing the numbers at
each time step for all species.

## See also

Other summary functions:
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
numbers <- getN(NS_sim)
numbers["1972", "Herring"]
#> [1] 1.49413e+11
# The above gave a huge number, because that included all the larvae.
# The number of Herrings between 10g and 1kg is much smaller.
numbers <- getN(NS_sim, min_w = 10, max_w = 1000)
numbers["1972", "Herring"]
#> [1] 4014916500
```
