# Calculate the total biomass of each species within a size range at each time step.

Calculates the total biomass through time within user defined size
limits. The default option is to use the size range starting at the size
specified by the `biomass_cutoff` species parameter, if it is set, or
else the full size range of each species. You can specify minimum and
maximum weight or length range for the species. Lengths take precedence
over weights (i.e. if both min_l and min_w are supplied, only min_l will
be used).

## Usage

``` r
getBiomass(object, use_cutoff = FALSE, ...)
```

## Arguments

- object:

  An object of class `MizerParams` or `MizerSim`.

- use_cutoff:

  If TRUE, the `biomass_cutoff` column in the species parameters is used
  as the minimum weight for each species (ignoring any size range
  arguments in `...`). If FALSE (default), the specified size range
  arguments are used, if provided, or the full size range of the species
  is used.

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

If called with a MizerParams object, a named vector with the biomass in
grams for each species in the model. If called with a MizerSim object,
an `ArrayTimeBySpecies` object (time x species) containing the biomass
in grams at each time step for all species.

## Details

When no size range arguments are provided, the function checks if the
`biomass_cutoff` column exists in the species parameters. If it does,
those values are used as the minimum weight for each species. For
species with NA values in `biomass_cutoff`, the default minimum weight
(smallest weight in the model) is used.

## See also

Other summary functions:
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
biomass <- getBiomass(NS_sim)
biomass["1972", "Herring"]
#> [1] 218218354800
biomass <- getBiomass(NS_sim, min_w = 10, max_w = 1000)
biomass["1972", "Herring"]
#> [1] 154290431041

# If species_params contains a `biomass_cutoff`` column, it can be used
# as the minimum weight when use_cutoff = TRUE
species_params(NS_sim@params)$biomass_cutoff <- 10
biomass <- getBiomass(NS_sim, use_cutoff = TRUE)  # Uses biomass_cutoff as min_w
biomass["1972", "Herring"]
#> [1] 154290431041
```
