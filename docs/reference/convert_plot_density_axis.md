# Express plotting data on the requested size axis

Converts the size coordinate of the plotting data to the requested axis
and, when the values are a density, multiplies them by the Jacobian that
restates them in the density measure that axis calls for (see
[`density_target_measure()`](https://sizespectrum.org/mizer/reference/density_target_measure.md)).
Values that are not a density are left alone.

## Usage

``` r
convert_plot_density_axis(
  plot_dat,
  params,
  size_axis,
  density_wrt = NA_character_,
  per_log_size = NULL,
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

- density_wrt:

  The measure the values are a density with respect to, see
  [density_measures](https://sizespectrum.org/mizer/reference/density_measures.md).
  `NA` (the default) means the values are not a density.

- per_log_size:

  Whether to express the values per logarithmic size. `NULL` (the
  default) keeps whichever the values already are.

- species_col:

  Name of the column identifying the species. Default is `"Species"`.

- value_col:

  Name or index of the value column. Defaults to the second column.

## Value

The plotting data with its size coordinate, and where called for its
values, expressed for the requested axis. The size coordinate is the
first column.

## Details

Anything involving a length is a per-species quantity, because the
weight-length relationship is, so rows whose species is not one of the
model's species — the "Total" row, for instance — are dropped when a
length is needed. Going from a density per size to one per logarithmic
size needs no length, and keeps those rows.
