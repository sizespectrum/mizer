# Rename gears

Changes the names of gears in a MizerParams object. This involves for
example changing the gear dimension names of selectivity and
catchability arrays appropriately.

## Usage

``` r
renameGear(params, replace, ...)
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

[`renameSpecies()`](https://sizespectrum.org/mizer/reference/renameSpecies.md)

## Examples

``` r
replace <- c(Industrial = "Trawl", Otter = "Beam_Trawl")
params <- renameGear(NS_params, replace)
gear_params(params)$gear
#>        Sprat, Trawl      Sandeel, Trawl       N.pout, Trawl    Herring, Pelagic 
#>             "Trawl"             "Trawl"             "Trawl"           "Pelagic" 
#>           Dab, Beam Whiting, Beam_Trawl          Sole, Beam Gurnard, Beam_Trawl 
#>              "Beam"        "Beam_Trawl"              "Beam"        "Beam_Trawl" 
#>        Plaice, Beam Haddock, Beam_Trawl     Cod, Beam_Trawl  Saithe, Beam_Trawl 
#>              "Beam"        "Beam_Trawl"        "Beam_Trawl"        "Beam_Trawl" 
```
