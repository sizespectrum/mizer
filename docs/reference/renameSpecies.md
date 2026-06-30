# Rename species

Changes the names of species in a MizerParams object. This involves for
example changing the species dimension names of rate arrays
appropriately.

## Usage

``` r
renameSpecies(params, replace, ...)
```

## Arguments

- params:

  A mizer params object

- replace:

  A named character vector, with new names as values, and old names as
  names.

- ...:

  Currently unused.

## Value

An object of type
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)

## See also

[`renameGear()`](https://sizespectrum.org/mizer/reference/renameGear.md)

## Examples

``` r
replace <- c(Cod = "Kabeljau", Haddock = "Schellfisch")
params <- renameSpecies(NS_params, replace)
species_params(params)$species
#>  [1] "Sprat"       "Sandeel"     "N.pout"      "Herring"     "Dab"        
#>  [6] "Whiting"     "Sole"        "Gurnard"     "Plaice"      "Schellfisch"
#> [11] "Kabeljau"    "Saithe"     
```
