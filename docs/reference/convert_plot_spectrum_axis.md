# Convert a weight-based spectrum to a length-based spectrum

A density with respect to weight is converted to a density with respect
to length with the Jacobian `dw/dl = b * w / l`. A density with respect
to logarithmic weight is instead converted with
`d log(w) / d log(l) = b`. This is the interface used by the
`power`-based spectrum plots; arrays carry their density measure
explicitly and use
[`convert_plot_density_axis()`](https://sizespectrum.org/mizer/reference/convert_plot_density_axis.md)
instead.

## Usage

``` r
convert_plot_spectrum_axis(
  plot_dat,
  params,
  size_axis,
  power,
  per_log_size = power == 2,
  species_col = "Species",
  value_col = 2
)
```

## Arguments

- plot_dat:

  A data frame of plotting data with a `w` column and a species column.

- params:

  A MizerParams object providing the weight-length parameters.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- power:

  The power of weight multiplying the number density.

- per_log_size:

  Whether the spectrum is a density with respect to logarithmic size
  rather than with respect to size. Defaults to `power == 2`, the only
  power for which this used to be the case.

- species_col:

  Name of the column identifying the species. Default is `"Species"`.

- value_col:

  Name or index of the value column. Defaults to the second column.

## Value

The plotting data with both its size coordinate and spectrum values
expressed for the requested axis.
