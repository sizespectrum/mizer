# The density measure a plot calls for

A density is expressed with respect to two independent choices: the size
variable, which follows `size_axis`, and whether it is per size or per
logarithmic size, which follows `per_log_size`. Plotting against a
length axis therefore turns a density with respect to weight into one
with respect to length, and a density with respect to logarithmic weight
into one with respect to logarithmic length.

## Usage

``` r
density_target_measure(density_wrt, size_axis, per_log_size = NULL)
```

## Arguments

- density_wrt:

  The measure the values are a density with respect to, see
  [density_measures](https://sizespectrum.org/mizer/reference/density_measures.md).

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- per_log_size:

  Whether to express the values per logarithmic size. `NULL` (the
  default) keeps whichever the values already are.

## Value

The density measure to express the values in, or `NA_character_` if the
values are not a density.
