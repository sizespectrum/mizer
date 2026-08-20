# The species parameter columns holding an observation

Internal helper for the calibration and matching functions. Observations
of type `to` live in the `<to>_observed` column of the species
parameters, alongside an optional `<to>_cutoff` column giving the
smallest size that the observation includes.

## Usage

``` r
observation_columns(to = c("biomass", "number"))
```

## Arguments

- to:

  The type of observation, either "biomass" or "number".

## Value

A list with entries `to`, `observed` and `cutoff`, the latter two giving
the names of the corresponding species parameter columns.
