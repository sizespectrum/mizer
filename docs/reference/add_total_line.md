# Add a total line to plotting data by summing over its series

The total has to be formed *after* the size coordinate has been
converted, not before. On a weight axis every series shares the model's
weight grid, so summing at equal weight and summing at equal position
are the same thing. On a length axis they are not: each species, and the
resource, converts weight to length with its own allometric
relationship, so at a given length the series sit at different weights
and their grids no longer coincide. The sum that means something there
is the sum at equal *length* — the number of organisms per unit length,
whatever they are — which is what this computes.

## Usage

``` r
add_total_line(
  plot_dat,
  x_var = names(plot_dat)[[1]],
  value_col = 2,
  by = NULL
)
```

## Arguments

- plot_dat:

  A data frame of plotting data with a size column, a value column and a
  `Species` column.

- x_var:

  Name of the size column. Defaults to the first column.

- value_col:

  Name or index of the value column. Defaults to the second.

- by:

  Names of further columns identifying separate plots, such as the time
  of an animation frame or the model of a comparison. A total is formed
  within each of their combinations.

## Value

`plot_dat` with the total appended as a series named `"Total"`.

## Details

Each series is interpolated onto the sorted union of all the size
coordinates, linearly in the logarithm of size, since the grid is
logarithmic. A series contributes nothing outside its own range. When
the series already share a grid — always on a weight axis, and on a
length axis whenever the weight-length parameters agree — the union is
that grid and the interpolation reproduces the values exactly, so
nothing is approximated in the cases where nothing needs to be.
