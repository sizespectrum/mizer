# Apply a `second_order_w` value to the current slot list

Internal helper that validates a `second_order_w` `value` (a single
logical, a single flux scheme name, or a named vector with entries
`flux` and/or `bin_average`) and returns the updated named list. Shared
by the `second_order_w<-` setter and by the model constructors (e.g.
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)),
which set the slot directly before the rest of the parameters are
computed so that the bin-averaged constructions pick up the flag,
without the setter's extra
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md)
call.

## Usage

``` r
validate_second_order_w(current, value)
```

## Arguments

- current:

  The current `second_order_w` slot list (with entries `flux` and
  `bin_average`).

- value:

  The value to apply, as described above.

## Value

The updated `second_order_w` list.
