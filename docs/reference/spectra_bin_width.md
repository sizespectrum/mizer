# Look up the size-bin widths for spectra data

Matches the weights in the plotting data to the model's full size grid
and returns the corresponding bin widths.

## Usage

``` r
spectra_bin_width(w, params)
```

## Arguments

- w:

  Numeric vector of weights.

- params:

  A MizerParams object.

## Value

A numeric vector of bin widths, one for each entry of `w`.
