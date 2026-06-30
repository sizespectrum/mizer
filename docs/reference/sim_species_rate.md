# Build a `MizerSim` rate getter that resolves the rate functions once

Like
[`sim_size_rate()`](https://sizespectrum.org/mizer/reference/sim_size_rate.md)
but for getters that return one value per species at each time step (a
time-by-species array), such as
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md) and
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md). By
default these use the initial effort, matching their `MizerParams`
counterparts.

## Usage

``` r
sim_species_rate(
  sim,
  time_range,
  target,
  slot,
  value_name,
  units = NULL,
  use_sim_effort = FALSE,
  ...
)
```

## Arguments

- sim:

  A `MizerSim` object.

- time_range:

  Passed to the sim iteration helper.

- target:

  Name of the rate to calculate (as in `params@rates_funcs`).

- slot:

  Name of the element to extract from the
  [`mizer_rates_subset()`](https://sizespectrum.org/mizer/reference/mizer_rates_subset.md)
  result (e.g. `"e_growth"` for the `EGrowth` rate).

- value_name, units:

  Metadata for the returned array.

- use_sim_effort:

  If `TRUE`, the saved effort at each time step is used; otherwise the
  initial effort is used (matching the behaviour of the corresponding
  `MizerParams` getter).

- ...:

  Passed on to the rate functions.

## Value

An `ArrayTimeBySpecies` object.
