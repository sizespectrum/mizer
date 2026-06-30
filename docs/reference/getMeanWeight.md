# Calculate the mean weight of the community

Calculates the mean weight of the community. This is simply the total
biomass of the community divided by the abundance in numbers. You can
specify minimum and maximum weight or length range for the species.
Lengths take precedence over weights (i.e. if both min_l and min_w are
supplied, only min_l will be used). You can also specify the species to
be used in the calculation.

## Usage

``` r
getMeanWeight(object, species = NULL, ...)
```

## Arguments

- object:

  A
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  or
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- species:

  The species to be selected. Optional. By default all target species
  are selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

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

A vector containing the mean weight of the community through time, or a
single value if called with a `MizerParams` object.

## See also

Other functions for calculating indicators:
[`getCommunitySlope()`](https://sizespectrum.org/mizer/reference/getCommunitySlope.md),
[`getMeanMaxWeight()`](https://sizespectrum.org/mizer/reference/getMeanMaxWeight.md),
[`getProportionOfLargeFish()`](https://sizespectrum.org/mizer/reference/getProportionOfLargeFish.md)

## Examples

``` r
mean_weight <- getMeanWeight(NS_sim)
years <- c("1967", "2010")
mean_weight[years]
#>      1967      2010 
#> 0.7959689 0.4835935 
getMeanWeight(NS_sim, species = c("Herring", "Sprat", "N.pout"))[years]
#>      1967      2010 
#> 0.6750526 0.5773987 
getMeanWeight(NS_sim, min_w = 10, max_w = 5000)[years]
#>     1967     2010 
#> 36.12859 50.54693 
getMeanWeight(NS_params)
#> [1] 1.174164
```
