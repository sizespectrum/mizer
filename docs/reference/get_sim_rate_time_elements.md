# Get selected saved time steps for a simulation rate

Internal helper used by `MizerSim` rate methods. If `time_range` is
missing, all saved simulation times are selected; otherwise the request
is delegated to
[`get_time_elements()`](https://sizespectrum.org/mizer/reference/get_time_elements.md).

## Usage

``` r
get_sim_rate_time_elements(sim, time_range)
```

## Arguments

- sim:

  A `MizerSim` object.

- time_range:

  A numeric or character vector of times.

## Value

A named logical vector indicating the selected saved time steps.
