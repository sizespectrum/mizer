# Y-axis limits an array's type calls for

Only a `"proportion"` has an opinion. A `"density"` is handled where the
size axis is converted, and a `"value"` needs nothing.

## Usage

``` r
array_ylim(x, ylim, log_y, values)
```

## Arguments

- x:

  A mizer array object.

- ylim:

  Numeric vector of length two, the limits the caller asked for.

- log_y:

  Whether the y axis is logarithmic.

- values:

  The values being plotted.

## Value

A numeric vector of length two.
