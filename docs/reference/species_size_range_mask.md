# Which size classes lie inside each species' own size range

A rate array is defined on the whole size grid, but a species only
occupies the part of it between its `w_min` and `w_max`. The values
outside that range are arithmetic rather than biology — the encounter
rate a 40 kg Sprat would have — so
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) and
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md) leave
them out unless asked for them with `all.sizes = TRUE`. This is the
single definition of which classes those are, so that the two cannot
disagree about what the array contains.

## Usage

``` r
species_size_range_mask(params, w, species)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object, or `NULL`.

- w:

  The size represented by each column, from
  [`get_ArraySpeciesBySize_w()`](https://sizespectrum.org/mizer/reference/get_ArraySpeciesBySize_w.md).
  Passed in rather than read from `params` because a bin-averaged array
  is drawn at the bin centres.

- species:

  Character vector of the species, one per row.

## Value

A logical matrix (species x size), `TRUE` where the size class lies
inside that species' size range.

## Details

A species that the model does not know about — a row whose name is not
in `species_params`, or any row at all when the array carries no
`params` — has no range to be outside of, so all of its size classes are
kept.
