# Reorder species

Reorders all core species axes and invokes registered extension hooks.
Every current species must occur exactly once in `species`.

## Usage

``` r
reorderSpecies(params, species, ...)
```

## Arguments

- params:

  A MizerParams object.

- species:

  Species names in the desired order, or a permutation of their current
  integer indices.

- ...:

  Unused.

## Value

The updated MizerParams object.

## See also

[`addSpecies()`](https://sizespectrum.org/mizer/reference/addSpecies.md),
[`removeSpecies()`](https://sizespectrum.org/mizer/reference/removeSpecies.md),
[`renameSpecies()`](https://sizespectrum.org/mizer/reference/renameSpecies.md)
