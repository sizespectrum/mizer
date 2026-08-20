# Observed yield of each species

The observed yield lives in the `yield_observed` column of the gear
parameter data frame, see
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
where it is given for each gear-species pair. This function adds the
observations up over the gears to give the total observed yield of each
species. With the `gear` argument you can restrict the sum to a subset
of the gears.

## Usage

``` r
get_yield_observed(params, gear = NULL)
```

## Arguments

- params:

  A MizerParams object

- gear:

  The gears whose observations are to be added up. Optional. By default
  all gears are included. A vector of gear names.

## Value

A numeric vector with one entry for each species, named by species,
holding the observed yield in grams per year, or `NA` for species
without an observation. `NULL` if the observations are not available: if
neither the gear parameters nor the species parameters have a
`yield_observed` column, or, when `gear` is given, if the gear
parameters have no such column.

## Details

Older models, and the examples in older versions of mizer, put
`yield_observed` into the species parameter data frame instead. That is
still accepted: a species that has no observation among the gear
parameters takes its value from the species parameters. Where both
tables give a value for a species, the gear parameters win. The species
parameter observation is a total over all gears, so it is ignored when
`gear` selects only some of the gears.
