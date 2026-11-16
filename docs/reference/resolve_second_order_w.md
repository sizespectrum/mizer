# Resolve a `second_order_w` value against the default scheme

Internal helper that validates a `second_order_w` `value` against the
default scheme of the defaults edition (see
[`default_second_order_w()`](https://sizespectrum.org/mizer/reference/default_second_order_w.md))
and returns the resulting named list. A `NULL` value gives the default
scheme, and a value that sets only one entry leaves the other at its
default. Used by the model constructors to work out the target `flux`
and `bin_average` entries before the rest of the model is built.

## Usage

``` r
resolve_second_order_w(value, edition)
```

## Arguments

- value:

  The value to resolve, as accepted by `second_order_w<-`, or `NULL`.

- edition:

  The defaults edition of the model that is being built.

## Value

The resolved `second_order_w` list.
