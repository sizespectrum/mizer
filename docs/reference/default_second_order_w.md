# Default numerical scheme of a defaults edition

Internal helper giving the `second_order_w` slot list that a new model
of the given
[`defaults_edition()`](https://sizespectrum.org/mizer/reference/defaults_edition.md)
gets when the user does not choose a scheme: the second-order scheme
(`flux = "van_leer"`, `bin_average = TRUE`) from edition 3 on and the
first-order scheme (`flux = "upwind"`, `bin_average = FALSE`) before.

## Usage

``` r
default_second_order_w(edition)
```

## Arguments

- edition:

  The defaults edition of the model.

## Value

The default `second_order_w` list.
