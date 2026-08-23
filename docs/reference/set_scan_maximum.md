# Record where each series attains its maximum

Sets the `at_max` and `max_value` attributes from the rows currently in
the object, so that subsetting cannot leave a stale maximum behind.

## Usage

``` r
set_scan_maximum(x)
```

## Arguments

- x:

  A MizerScan object.

## Value

`x` with the two attributes set.
