# Designate species as background species

Marks the specified set of species as background species by setting the
`is_background` column in their species parameters to `TRUE`. Background
species are handled differently in plots (displayed in grey) and their
abundances can be automatically adjusted to keep the community close to
the Sheldon spectrum (see `adjustBackgroundSpecies()` in the
mizerExperimental package).

## Usage

``` r
markBackground(object, species = NULL)
```

## Arguments

- object:

  An object of class
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  or
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md).

- species:

  The species to be selected. Optional. By default all target species
  are selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

## Value

An object of the same class as the `object` argument

## See also

[`removeBackgroundSpecies()`](https://sizespectrum.org/mizer/reference/removeBackgroundSpecies.md)

## Examples

``` r
params <- markBackground(NS_params,
                         species = c("Sprat", "Sandeel", "N.pout"))
any(species_params(params)$is_background)
#> [1] TRUE
```
