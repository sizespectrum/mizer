# Calculate the slope of the community abundance

Calculates the slope of the community abundance by performing a linear
regression on the logged total numerical abundance at weight and logged
weights (natural logs, not log to base 10, are used). You can specify
minimum and maximum weight or length range for the species. Lengths take
precedence over weights (i.e. if both min_l and min_w are supplied, only
min_l will be used). You can also specify the species to be used in the
calculation.

## Usage

``` r
getCommunitySlope(object, species = NULL, biomass = TRUE, ...)
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

- biomass:

  Boolean. If TRUE (default), the abundance is based on biomass, if
  FALSE the abundance is based on numbers.

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

A data.frame with columns slope, intercept and the coefficient of
determination R^2 (and a time step column when called with a `MizerSim`
object).

## See also

Other functions for calculating indicators:
[`getMeanMaxWeight()`](https://sizespectrum.org/mizer/reference/getMeanMaxWeight.md),
[`getMeanWeight()`](https://sizespectrum.org/mizer/reference/getMeanWeight.md),
[`getProportionOfLargeFish()`](https://sizespectrum.org/mizer/reference/getProportionOfLargeFish.md)

## Examples

``` r
# Slope based on biomass, using all species and sizes
slope_biomass <- getCommunitySlope(NS_sim)
slope_biomass[1, ] # in 1976
#>          slope intercept        r2
#> 1967 -0.830658  25.43745 0.8878355
slope_biomass[idxFinalT(NS_sim), ] # in 2010
#>           slope intercept       r2
#> 2010 -0.7925937  25.61399 0.935702

# Slope based on numbers, using all species and sizes
slope_numbers <- getCommunitySlope(NS_sim, biomass = FALSE)
slope_numbers[1, ] # in 1976
#>          slope intercept        r2
#> 1967 -1.830658  25.43745 0.9746487

# Slope based on biomass, using all species and sizes between 10g and 1000g
slope_biomass <- getCommunitySlope(NS_sim, min_w = 10, max_w = 1000)
slope_biomass[1, ] # in 1976
#>          slope intercept        r2
#> 1967 -1.526323  30.00691 0.9337066

# Slope based on biomass, using only demersal species and
# sizes between 10g and 1000g
dem_species <- c("Dab","Whiting", "Sole", "Gurnard", "Plaice",
                 "Haddock", "Cod", "Saithe")
slope_biomass <- getCommunitySlope(NS_sim, species = dem_species,
                                   min_w = 10, max_w = 1000)
slope_biomass[1, ] # in 1976
#>           slope intercept        r2
#> 1967 -0.9704957  26.69254 0.8095734

getCommunitySlope(NS_params)
#>       slope intercept        r2
#> 1 -0.782225  25.40779 0.8722251
```
