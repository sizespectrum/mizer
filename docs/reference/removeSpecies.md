# Remove species

This function simply removes all entries from the MizerParams object
that refer to the selected species. It does not recalculate the steady
state for the remaining species or retune their reproductive efficiency.

## Usage

``` r
removeSpecies(params, species, ...)
```

## Arguments

- params:

  A mizer params object for the original system.

- species:

  The species to be removed. A vector of species names, or a numeric
  vector of species indices, or a logical vector indicating for each
  species whether it is to be removed (TRUE) or not.

- ...:

  Currently unused.

## Value

An object of type
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)

## Details

If a gear was targeting only the removed species, then this function
will NOT remove that gear. If you want to also remove that gear then you
can do that by calling
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md).

## See also

[`addSpecies()`](https://sizespectrum.org/mizer/reference/addSpecies.md),
[`renameSpecies()`](https://sizespectrum.org/mizer/reference/renameSpecies.md)

## Examples

``` r
params <- NS_params
species_params(params)$species
#>  [1] "Sprat"   "Sandeel" "N.pout"  "Herring" "Dab"     "Whiting" "Sole"   
#>  [8] "Gurnard" "Plaice"  "Haddock" "Cod"     "Saithe" 
params <- removeSpecies(params, c("Cod", "Haddock"))
species_params(params)$species
#>  [1] "Sprat"   "Sandeel" "N.pout"  "Herring" "Dab"     "Whiting" "Sole"   
#>  [8] "Gurnard" "Plaice"  "Saithe" 
```
