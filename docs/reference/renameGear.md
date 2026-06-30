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
#>  [1] "Trawl"      "Trawl"      "Trawl"      "Pelagic"    "Beam"      
#>  [6] "Beam_Trawl" "Beam"       "Beam_Trawl" "Beam"       "Beam_Trawl"
#> [11] "Beam_Trawl" "Beam_Trawl"
```
