# Adjust the size grid

This function adjusts the size grid in a
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.md)
object to the desired minimum and maximum size. It can both expand and
truncate the grid. If the grid is truncated, any data outside the new
grid is discarded. A warning is issued if there is non-negligible
biomass in the discarded size bins.

## Usage

``` r
adjustSizeGrid(params, ...)

# S3 method for class 'MizerParams'
adjustSizeGrid(
  params,
  new_min_w = min(params@species_params$w_min),
  new_max_w = max(params@species_params$w_max),
  new_min_w_pp = min(params@w_full),
  preserve_species = params@species_params$species,
  tol = 1e-06,
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.md)
  object.

- ...:

  Additional arguments.

- new_min_w:

  The new minimum size in the grid. Defaults to the minimum egg size of
  all species.

- new_max_w:

  The new maximum size in the grid. Defaults to the maximum asymptotic
  size of all species.

- new_min_w_pp:

  The new minimum size of the resource spectrum. Defaults to the current
  minimum of `w_full`.

- preserve_species:

  A vector of species names for which all rate arrays should be copied
  over to the new params object rather than being re-calculated from the
  species parameters. If missing, all species are preserved.

- tol:

  A numeric value specifying the tolerance for truncation losses. The
  following checks are made separately for each species and a warning is
  raised, listing the affected species, if the lost fraction exceeds
  this value for any of them: the fraction of the species' biomass lost,
  the fraction of the diet of the smallest individuals of the species
  lost to resource truncation, and the fraction of the diet of the
  largest individuals of the species lost to resource truncation.
  Defaults to `1e-6`.

## Value

A new
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.md)
object with the updated size grid.
