# Apply a species rate function over saved simulation times

Internal helper used by `MizerSim` rate methods whose one-time result is
a named vector with one value for each species.

## Usage

``` r
get_species_time_rate_from_sim(
  sim,
  time_range,
  rate_fun,
  value_name,
  units = NULL
)
```

## Arguments

- sim:

  A `MizerSim` object.

- time_range:

  A numeric or character vector of times.

- rate_fun:

  A function accepting a single simulation slice as returned by
  [`get_sim_rate_slice()`](https://sizespectrum.org/mizer/reference/get_sim_rate_slice.md).

- value_name:

  Name of the value stored in the returned array.

- units:

  Optional units of the value stored in the returned array.

## Value

An `ArrayTimeBySpecies` object with dimensions time x species.
