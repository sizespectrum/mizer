# Expand the size grid

**\[deprecated\]** This function expands the size grid in a
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.md)
object to the desired min and max size, preserving all existing species.
The function is deprecated because you can achieve the same more
flexibly with
[`adjustSizeGrid()`](https://sizespectrum.org/mizer/reference/adjustSizeGrid.md).

## Usage

``` r
expandSizeGrid(params, ...)

# S3 method for class 'MizerParams'
expandSizeGrid(
  params,
  new_min_w = min(params@w),
  new_max_w = max(params@w),
  preserve_species = params@species_params$species,
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.md)
  object.

- ...:

  Additional arguments (currently unused).

- new_min_w:

  The new minimum size in the grid. Defaults to the current minimum.

- new_max_w:

  The new maximum size in the grid. Defaults to the current maximum.

- preserve_species:

  A vector of species names for which all rate arrays should be copied
  over to the new params object rather than being re-calculated from the
  species parameters. If missing, all species are preserved.

## Value

A new
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.md)
object with the updated size grid.
