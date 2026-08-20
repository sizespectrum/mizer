# Convert plotting data from weight to length

When `size_axis = "l"`, adds a length column `l` computed from the
weight column `w` using the weight-length relationship of each line, see
[`plot_length_params()`](https://sizespectrum.org/mizer/reference/plot_length_params.md).
Rows with no such relationship are dropped. For `size_axis = "w"` the
data is returned unchanged.

## Usage

``` r
convert_plot_size_axis(
  plot_dat,
  params,
  size_axis,
  species_col = "Species",
  drop_w = TRUE
)
```

## Arguments

- plot_dat:

  A data frame of plotting data with a `w` column and a species column.

- params:

  A MizerParams object providing the weight-length parameters.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- species_col:

  Name of the column identifying the species. Default is `"Species"`.

- drop_w:

  If `TRUE` (default), the `w` column is dropped once `l` has been
  computed.

## Value

The plotting data, with a length column `l` added (and `w` optionally
dropped) when `size_axis = "l"`.
