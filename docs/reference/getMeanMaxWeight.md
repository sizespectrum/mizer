# Calculate the mean maximum weight of the community

Calculates the mean maximum weight of the community. This can be
calculated by numbers or biomass. The calculation is the sum of the
`w_inf` \* abundance of each species, divided by the total abundance
community, where abundance is either in biomass or numbers. You can
specify minimum and maximum weight or length range for the species.
Lengths take precedence over weights (i.e. if both min_l and min_w are
supplied, only min_l will be used). You can also specify the species to
be used in the calculation.

## Usage

``` r
getMeanMaxWeight(object, species = NULL, measure = "both", ...)
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

- measure:

  The measure to return. Can be 'numbers', 'biomass' or 'both'

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

Depends on the `measure` argument. If `measure = “both”` then you get a
matrix with two columns, one with values by numbers, the other with
values by biomass at each saved time step (or a named vector with two
entries for `MizerParams`). If `measure = “numbers”` or `“biomass”` you
get a vector of the respective values at each saved time step (or a
single value for `MizerParams`).

## See also

Other functions for calculating indicators:
[`getCommunitySlope()`](https://sizespectrum.org/mizer/reference/getCommunitySlope.md),
[`getMeanWeight()`](https://sizespectrum.org/mizer/reference/getMeanWeight.md),
[`getProportionOfLargeFish()`](https://sizespectrum.org/mizer/reference/getProportionOfLargeFish.md)

## Examples

``` r
mmw <- getMeanMaxWeight(NS_sim)
years <- c("1967", "2010")
mmw[years, ]
#>      mmw_numbers mmw_biomass
#> 1967    2640.075   10980.365
#> 2010    2716.657    7940.787
getMeanMaxWeight(NS_sim, species=c("Herring","Sprat","N.pout"))[years, ]
#>      mmw_numbers mmw_biomass
#> 1967    121.2163    226.2291
#> 2010    122.7208    222.6168
getMeanMaxWeight(NS_sim, min_w = 10, max_w = 5000)[years, ]
#>      mmw_numbers mmw_biomass
#> 1967    1313.207    3444.342
#> 2010    2316.711    6742.697
getMeanMaxWeight(NS_params)
#> mmw_numbers mmw_biomass 
#>    2335.168    4583.528 
```
