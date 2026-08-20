# Match a quantity to observations species by species

Internal implementation shared by
[`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md)
and
[`matchNumbers()`](https://sizespectrum.org/mizer/reference/matchNumbers.md).
Multiplies the abundance density of each selected species at all sizes
by the factor that brings the modelled quantity onto the observation.
Species that were not selected, or that have no positive observation,
are left alone.

## Usage

``` r
match_to(params, species = NULL, to = c("biomass", "number"), fname)
```

## Arguments

- params:

  A MizerParams object.

- species:

  The species to be affected, in any of the forms accepted by
  [`valid_species_arg()`](https://sizespectrum.org/mizer/reference/valid_species_arg.md).

- to:

  The type of observation, either "biomass" or "number".

- fname:

  The name of the calling function, used when reporting that the model
  has been moved off its steady state.

## Value

A MizerParams object.
