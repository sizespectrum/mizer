# Interpolate a series linearly in the logarithm of size

The size grid is logarithmic, so the interpolation is too. A series of a
single point can only speak for the coordinate it sits at, and a
coordinate that is not positive has no logarithm, so those two cases
fall back to matching the coordinates exactly and to a linear
interpolation respectively.

## Usage

``` r
interpolate_in_log_size(x, y, xout)
```

## Arguments

- x, y:

  The coordinates and values of the series.

- xout:

  The coordinates to interpolate onto.

## Value

The interpolated values, `NA` where `xout` is outside the range of `x`.
