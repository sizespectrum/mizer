# Get available energy

**\[deprecated\]**

This is deprecated and is no longer used by the mizer project() method.
Calculates the amount \\E\_{a,i}(w)\\ of food exposed to each predator
as a function of predator size.

## Usage

``` r
getPhiPrey(object, n, n_pp, ...)
```

## Arguments

- object:

  An
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- n:

  A matrix of species abundances (species x size)

- n_pp:

  A vector of the background abundance by size

- ...:

  Other arguments (currently unused)

## Value

A two dimensional array (predator species x predator size) equal to
`getEncounter(object, n, n_pp) / getSearchVolume(object)`.

## See also

[`project()`](https://sizespectrum.org/mizer/reference/project.md)
