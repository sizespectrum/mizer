# Convert plotting data from weight to length

When `size_axis = "l"`, adds a length column `l` computed from the
weight column `w` using each species' weight-length relationship. For
`size_axis = "w"` the data is returned unchanged.

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
