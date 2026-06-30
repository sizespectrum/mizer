# Calculate the proportion of large fish

Calculates the proportion of large fish in a `MizerSim` or `MizerParams`
object within user defined size limits. The default option is to use the
whole size range. You can specify minimum and maximum size ranges for
the species and also the threshold size for large fish. Sizes can be
expressed as weight or length. Lengths take precedence over weights
(i.e. if both `min_l` and `min_w` are supplied, only `min_l` will be
used, and if `threshold_l` is supplied it takes precedence over
`threshold_w`). You can also specify the species to be used in the
calculation. This function can be used to calculate the Large Fish
Index. The proportion is based on either abundance or biomass.

## Usage

``` r
getProportionOfLargeFish(
  object,
  species = NULL,
  threshold_w = 100,
  threshold_l = NULL,
  biomass_proportion = TRUE,
  ...
)
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

- threshold_w:

  The weight used as the cutoff between large and small fish. Default
  value is 100.

- threshold_l:

  The length used as the cutoff between large and small fish. If
  supplied, this takes precedence over `threshold_w`.

- biomass_proportion:

  A boolean value. If TRUE the proportion calculated is based on
  biomass, if FALSE it is based on numbers of individuals. Default is
  TRUE.

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

A vector containing the proportion of large fish through time, or a
single value if called with a `MizerParams` object.

## See also

Other functions for calculating indicators:
[`getCommunitySlope()`](https://sizespectrum.org/mizer/reference/getCommunitySlope.md),
[`getMeanMaxWeight()`](https://sizespectrum.org/mizer/reference/getMeanMaxWeight.md),
[`getMeanWeight()`](https://sizespectrum.org/mizer/reference/getMeanWeight.md)

## Examples

``` r
lfi <- getProportionOfLargeFish(NS_sim, min_w = 10, max_w = 5000,
                                threshold_w = 500)
years <- c("1972", "2010")
lfi[years]
#>      1972      2010 
#> 0.1425938 0.2573979 
getProportionOfLargeFish(NS_sim)[years]
#>      1972      2010 
#> 0.3115163 0.4441647 
getProportionOfLargeFish(NS_sim, species=c("Herring","Sprat","N.pout"))[years]
#>       1972       2010 
#> 0.06371869 0.21844009 
getProportionOfLargeFish(NS_sim, min_w = 10, max_w = 5000)[years]
#>      1972      2010 
#> 0.3474937 0.5662073 
getProportionOfLargeFish(NS_sim, min_w = 10, max_w = 5000,
    threshold_w = 500, biomass_proportion = FALSE)[years]
#>       1972       2010 
#> 0.00441820 0.01212436 
getProportionOfLargeFish(NS_params)
#> [1] 2.312485e-07
```
