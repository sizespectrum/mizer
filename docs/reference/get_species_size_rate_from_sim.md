# Apply a species-by-size rate function over saved simulation times

Internal helper used by `MizerSim` rate methods whose one-time result is
an `ArraySpeciesBySize`. The helper applies the supplied rate function
to each selected time slice, stacks the results, and restores the
appropriate mizer array class when dimensions have not been dropped.

## Usage

``` r
get_species_size_rate_from_sim(
  sim,
  time_range,
  drop,
  rate_fun,
  value_name,
  units = NULL,
  representation = "point"
)
```

## Arguments

- sim:

  A `MizerSim` object.

- time_range:

  A numeric or character vector of times.

- drop:

  If `TRUE`, dimensions of length 1 are dropped from the result.

- rate_fun:

  A function accepting a single simulation slice as returned by
  [`get_sim_rate_slice()`](https://sizespectrum.org/mizer/reference/get_sim_rate_slice.md).

- value_name:

  Name of the value stored in the returned array.

- units:

  Optional units of the value stored in the returned array.

## Value

A time x species x size array, possibly with dimensions dropped.
