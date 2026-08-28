# Calculate the mean size of the community

`getMeanWeight()` calculates the mean weight of the community. This is
simply the total biomass of the community divided by the abundance in
numbers. `getMeanLength()` calculates the mean length of the community,
i.e. the total length of all individuals divided by their number, where
the length of an individual of weight \\w\\ of species \\i\\ is obtained
from the length-weight relationship \\w = a_i l^{b_i}\\ with the species
parameters `a` and `b`.

## Usage

``` r
getMeanWeight(object, species = NULL, ...)

getMeanLength(object, species = NULL, ...)
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

A vector containing the mean weight (in grams) or the mean length (in
cm) of the community through time, or a single value if called with a
`MizerParams` object.

## Details

The length-weight parameters `a` and `b` are taken from the species
parameter data frame, where they are given the defaults `a = 0.01` and
`b = 3` when a model is created. `getMeanLength()` gives an error if a
model is so old that it has no such columns.

You can specify minimum and maximum weight or length for the included
size range. Lengths take precedence over weights (i.e. if both min_l and
min_w are supplied, only min_l will be used). You can also specify the
species to be used in the calculation.

You will usually want to give a minimum size. Over the full size range
the community is dominated in numbers by the smallest individuals, so
both means describe the larvae rather than the fish, and in the case of
`getMeanLength()` they do so with a length-weight relationship that was
fitted to observed fish and is being extrapolated far below the sizes it
describes. A minimum size also makes the indicator comparable to one
calculated from survey data, which sees only the sizes the gear catches.

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

mean_length <- getMeanLength(NS_sim)
mean_length[years]
#>     1967     2010 
#> 1.822872 1.625172 
# Only a couple of centimetres, because the larvae outnumber everything
# else. Give a minimum size to get an indicator about fish:
getMeanLength(NS_sim, min_l = 10)[years]
#>     1967     2010 
#> 17.50719 15.77286 
getMeanLength(NS_sim, min_w = 10, max_w = 5000)[years]
#>     1967     2010 
#> 19.70678 18.40549 
getMeanLength(NS_sim, species = c("Herring", "Sprat", "N.pout"),
              min_l = 10)[years]
#>     1967     2010 
#> 15.44524 14.48621 
getMeanLength(NS_sim@params, min_l = 10)
#> [1] 17.50719
```
