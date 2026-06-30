# Get array indices for a time range in a MizerSim object

Internal helper to select the saved time points whose times lie between
the smallest and largest values in `time_range`, inclusive.

## Usage

``` r
get_time_elements(sim, time_range, slot_name = "n")
```

## Arguments

- sim:

  A MizerSim object.

- time_range:

  A numeric or character vector of times. Only the range of values
  matters, so all saved times between `min(time_range)` and
  `max(time_range)` are selected.

- slot_name:

  Obsolete, kept only for backward compatibility with early versions
  where different time-based slots could have different time grids.
  Leave at the default.

## Value

A named logical vector, with one entry for each saved time in `sim`,
indicating whether that time lies in the requested range.
