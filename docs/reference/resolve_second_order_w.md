# Resolve a `second_order_w` value against the default scheme

Internal helper that validates a `second_order_w` `value` against the
default first-order slot (`flux = "upwind"`, `bin_average = FALSE`) and
returns the resulting named list. Used by the model constructors to work
out the target `flux` and `bin_average` entries before the rest of the
model is built.

## Usage

``` r
resolve_second_order_w(value)
```

## Arguments

- value:

  The value to resolve, as accepted by `second_order_w<-`.

## Value

The resolved `second_order_w` list.
