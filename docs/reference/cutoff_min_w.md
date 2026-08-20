# The minimum weights given by a cutoff species parameter

Internal helper for
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
and for the calibration and matching functions. Returns the
`<to>_cutoff` column of the species parameters, with any NAs replaced by
the smallest weight in the model. If the model has no such column at
all, the smallest weight is used for every species, so that the whole
size range is counted.

## Usage

``` r
cutoff_min_w(params, to = c("biomass", "number"))
```

## Arguments

- params:

  A MizerParams object.

- to:

  The type of observation, either "biomass" or "number".

## Value

A numeric vector with one minimum weight for each species.
