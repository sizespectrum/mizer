# Set species interaction matrix

Set species interaction matrix

## Usage

``` r
setInteraction(params, interaction = NULL, ...)

interaction_matrix(params)

interaction_matrix(params) <- value
```

## Arguments

- params:

  MizerParams object

- interaction:

  Optional interaction matrix of the species (predator species x prey
  species). By default all entries are 1. See "Setting interaction
  matrix" section below.

- ...:

  Unused

- value:

  An interaction matrix

## Value

`setInteraction`: A MizerParams object with updated interaction matrix

`interaction_matrix()`: The interaction matrix (predator species x prey
species)

## Setting interaction matrix

You do not need to specify an interaction matrix. If you do not, then
the predator-prey interactions are purely determined by the size of
predator and prey and totally independent of the species of predator and
prey.

The interaction matrix \\\theta\_{ij}\\ modifies the interaction of each
pair of species in the model. This can be used for example to allow for
different spatial overlap among the species. The values in the
interaction matrix are used to scale the encountered food and predation
mortality (see on the website [the section on predator-prey encounter
rate](https://sizespectrum.org/mizer/articles/model_description.html#sec:pref)
and on [predation
mortality](https://sizespectrum.org/mizer/articles/model_description.html#mortality)).
The first index refers to the predator species and the second to the
prey species.

The interaction matrix is used when calculating the food encounter rate
in
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
and the predation mortality rate in
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md).
Its entries are dimensionless numbers. If all the values in the
interaction matrix are equal then predator-prey interactions are
determined entirely by size-preference.

This function checks that the supplied interaction matrix is valid and
then stores it in the `interaction` slot of the `params` object.

The order of the columns and rows of the `interaction` argument should
be the same as the order in the species params data frame in the
`params` object. If you supply a named array then the function will
check the order and message if it is different before ignoring the
supplied dimnames. If you supply only column names then these are also
used as the row names. One way of creating your own interaction matrix
is to enter the data using a spreadsheet program and saving it as a .csv
file. The data can then be read into R using the command
[`read.csv()`](https://rdrr.io/r/utils/read.table.html).

The interaction of the species with the resource are set via a column
`interaction_resource` in the `species_params` data frame. By default
this column is set to all 1s.

## See also

Other functions for setting parameters:
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
[`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md),
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)

## Examples

``` r
params <- newTraitParams(no_sp = 3)
inter <- getInteraction(params)
#> Warning: `getInteraction()` was deprecated in mizer 2.4.0.
#> ℹ Please use `interaction_matrix()` instead.
inter[1, 2:3] <- 0
params <- setInteraction(params, interaction = inter)
getInteraction(params)
#>         prey
#> predator 1 2 3
#>        1 1 0 0
#>        2 1 1 1
#>        3 1 1 1
```
