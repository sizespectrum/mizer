# Remove all background species

Removes all species that have been marked as background species with
[`markBackground()`](https://sizespectrum.org/mizer/reference/markBackground.md).

## Usage

``` r
removeBackgroundSpecies(params)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

## Value

A
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
object with background species removed

## Details

This is just a shorthand for
`removeSpecies(params, species_params(params)$is_background)`

## See also

[`markBackground()`](https://sizespectrum.org/mizer/reference/markBackground.md)

## Examples

``` r
params <- markBackground(NS_params,
                         species = c("Sprat", "Sandeel", "N.pout"))
params <- removeBackgroundSpecies(params)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
species_params(params)$species
#> [1] "Herring" "Dab"     "Whiting" "Sole"    "Gurnard" "Plaice"  "Haddock"
#> [8] "Cod"     "Saithe" 
```
