# Get diet of predator at size, resolved by prey species

Calculates the rate at which a predator of a particular species and size
consumes biomass of each prey species, resource, and other components of
the ecosystem. Returns either the rates in grams/year or the proportion
of the total consumption rate.

## Usage

``` r
getDiet(object, proportion = TRUE, ...)
```

## Arguments

- object:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  or
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  object.

- proportion:

  If TRUE (default) the function returns the diet as a proportion of the
  total consumption rate. If FALSE it returns the consumption rate in
  grams per year.

- ...:

  Additional arguments that depend on the class of `object`.

  **For a
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object:**

  `n`

  :   A matrix of species abundances (species x size). Defaults to the
      initial abundances stored in `object`.

  `n_pp`

  :   A vector of the resource abundance by size. Defaults to the
      initial resource abundance stored in `object`.

  `n_other`

  :   A named list of the abundances of other dynamical components.
      Defaults to the initial values stored in `object`.

  **For a
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  object:**

  `time_range`

  :   The time range over which to return the diet. Either a vector of
      values, a vector of min and max time, or a single value. Defaults
      to the whole time range of the simulation.

  `drop`

  :   If `TRUE` then any dimension of length 1 is removed from the
      returned array.

## Value

- `MizerParams`: An array (predator species x predator size x (prey
  species + resource + other components)). Dimnames are "predator", "w",
  and "prey".

- `MizerSim`: A four-dimensional array (time x predator species x
  predator size x prey) with the diet at each selected saved time step.
  If `drop = TRUE` then dimensions of length 1 are removed.

## Details

The rates \\D\_{ij}(w)\\ at which a predator of species \\i\\ and size
\\w\\ consumes biomass from prey species \\j\\ are calculated from the
predation kernel \\\phi_i(w, w_p)\\, the search volume \\\gamma_i(w)\\,
the feeding level \\f_i(w)\\, the species interaction matrix
\\\theta\_{ij}\\ and the prey abundance density \\N_j(w_p)\\: \$\$
D\_{ij}(w, w_p) = (1-f_i(w)) \gamma_i(w) \theta\_{ij} \int N_j(w_p)
\phi_i(w, w_p) w_p dw_p. \$\$ The prey index \\j\\ runs over all species
and the resource.

Extra columns are added for the external encounter rate and for any
extra ecosystem components in your model for which you have defined an
encounter rate function. These encounter rates are multiplied by
\\1-f_i(w)\\ to give the rate of consumption of biomass from these extra
components.

This function performs the same integration as
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
but does not aggregate over prey species, and multiplies by \\1-f_i(w)\\
to get the consumed biomass rather than the available biomass. Outside
the range of sizes for a predator species the returned rate is zero.

## See also

[`plotDiet()`](https://sizespectrum.org/mizer/reference/plotDiet.md)

Other summary functions:
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
diet <- getDiet(NS_params)
str(diet)
#>  num [1:12, 1:100, 1:14] 8.94e-18 6.86e-19 3.46e-18 1.75e-09 1.12e-17 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ predator: chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w       : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>   ..$ prey    : chr [1:14] "Sprat" "Sandeel" "N.pout" "Herring" ...
# \donttest{
# For a MizerSim the diet is returned at each saved time step
sim <- project(NS_params, t_max = 20, effort = 0.5)
# Diet at the saved time steps over years 15 - 20
diet <- getDiet(sim, time_range = c(15, 20))
str(diet)
#>  num [1:6, 1:12, 1:100, 1:14] 8.42e-18 1.37e-17 8.76e-18 7.98e-18 7.49e-18 ...
#>  - attr(*, "dimnames")=List of 4
#>   ..$ time    : chr [1:6] "15" "16" "17" "18" ...
#>   ..$ predator: chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w       : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>   ..$ prey    : chr [1:14] "Sprat" "Sandeel" "N.pout" "Herring" ...
# }
```
