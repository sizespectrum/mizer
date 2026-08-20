# Check that `per_log_size` applies to a mizer array

Expressing values per logarithmic size only means anything for a
density, so asking for it on anything else is an argument error rather
than something to be quietly ignored — which is what `...` used to do
with it.

## Usage

``` r
check_per_log_size(x, per_log_size)
```

## Arguments

- x:

  A mizer array object.

- per_log_size:

  The `per_log_size` argument of the plot method.

## Value

`per_log_size`, invisibly, if it applies.
