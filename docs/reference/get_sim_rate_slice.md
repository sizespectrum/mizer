# Extract one saved simulation state for a rate calculation

Internal helper used by `MizerSim` rate methods to rebuild the
single-time inputs expected by `MizerParams` rate methods.

## Usage

``` r
get_sim_rate_slice(sim, time_idx)
```

## Arguments

- sim:

  A `MizerSim` object.

- time_idx:

  Integer index of the saved time step to extract.

## Value

A list with entries `n`, `n_pp`, `n_other`, `effort`, and `t`.
