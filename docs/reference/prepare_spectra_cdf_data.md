# Integrate spectra data into a cumulative distribution

Multiplies the spectra density by the size-bin widths and forms the
cumulative sum over size for each species, optionally normalising each
curve to end at 1.

## Usage

``` r
prepare_spectra_cdf_data(plot_dat, params, normalise = TRUE)
```

## Arguments

- plot_dat:

  Spectra plotting data with a `w` column, a `Species` column and the
  value in the second column.

- params:

  A MizerParams object, used for the size-bin widths.

- normalise:

  If `TRUE` (default), divide each species' curve by its maximum.

## Value

The plotting data with the value column replaced by its cumulative
distribution over size.
