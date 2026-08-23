# Put two series on a common size grid before comparing them

A relative difference can only be formed where both series have a value.
On a weight axis they always do: both models share the size grid, so the
two frames match up row for row. On a length axis they need not, because
each model converts weight to length with its own allometric
relationship, so the same weight grid lands on different lengths.
Matching the frames by equality of the size coordinate then throws away
nearly every point — an inner join on two grids that merely overlap
keeps only their exact coincidences.

## Usage

``` r
interpolate_relative_frames(frame1, frame2, x_var, y_var, by_vars)
```

## Arguments

- frame1, frame2:

  Data frames of prepared plotting data, sharing their variable names.

- x_var:

  Name of the size column.

- y_var:

  Name of the value column.

- by_vars:

  Names of the columns identifying a series, typically the species and
  the legend group.

## Value

A data frame with the `by_vars`, the size column, and the two value
columns named `<y_var>.x` and `<y_var>.y`, holding only the series
present in both frames.

## Details

Each series is therefore interpolated, linearly in the logarithm of size
since the grid is logarithmic, onto the sorted union of the two sets of
coordinates, restricted to the interval both series cover. Outside that
interval one of them would have to be extrapolated, which is not a
comparison but a guess. When the two grids already coincide the union is
that grid, the overlap is all of it, and the interpolation reproduces
the values exactly, so the matching case is unchanged.
