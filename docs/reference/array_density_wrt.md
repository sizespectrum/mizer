# The density measure of a mizer array

The bridge from the array metadata into the density machinery of the
plots. Mizer arrays are indexed by the model's weight grid, so a stored
density is always a density with respect to weight; the other measures
in
[density_measures](https://sizespectrum.org/mizer/reference/density_measures.md)
arise only for quantities that the spectrum plots compute on the fly,
such as a density per logarithmic weight.

## Usage

``` r
array_density_wrt(x)
```

## Arguments

- x:

  A mizer array object.

## Value

`"w"` if the array holds a density, otherwise `NA_character_`.
