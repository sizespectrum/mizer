# Build a `MizerSim` rate getter that resolves the rate functions once

Internal helper capturing the pattern shared by the `MizerSim` rate
getters that return a species-by-size array. It validates the params and
resolves the rate functions a single time, then for each saved time step
calculates only the required `target` rate with
[`mizer_rates_subset()`](https://sizespectrum.org/mizer/reference/mizer_rates_subset.md)
and extracts the element named `slot` from the result.

## Usage

``` r
sim_size_rate(
  sim,
  time_range,
  drop,
  target,
  slot,
  value_name,
  units = NULL,
  use_sim_effort = FALSE,
  representation = "point",
  ...
)
```

## Arguments

- sim:

  A `MizerSim` object.

- time_range:

  Passed to the sim iteration helper.

- drop:

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

An `ArrayTimeBySpeciesBySize` object (or a reduced array if `drop`).
