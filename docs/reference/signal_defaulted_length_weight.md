# Report a selectivity built from a defaulted weight-length relationship

Internal helper for
[`calc_selectivity()`](https://sizespectrum.org/mizer/reference/calc_selectivity.md).
A length-based selectivity function converts the lengths in
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
to weights with the species' `a` and `b`. When those were filled in by
mizer rather than supplied, the selectivity curve sits at weights mizer
invented, and the fishing mortality is wrong in a way nothing else
reveals. So it is reported where the conversion happens and not only
where the default is filled in.

## Usage

``` r
signal_defaulted_length_weight(params, species)
```

## Arguments

- params:

  A MizerParams object.

- species:

  The species whose selectivity is set from their length.

## Value

`NULL`, invisibly. Called for the report it raises.
