# Get encounter rate of predator at size, resolved by food source

Calculates the rate at which a predator of a particular species and size
encounters biomass of each prey species, of the resource, and of the
other food sources in the model. Summed over food sources this is the
encounter rate returned by
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md).
Multiplied by \\1-f_i(w)\\ it is the diet returned by
`getDiet(proportion = FALSE)`, and that is what this function is mainly
for: it is how
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md) and
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md)
split the encounter rate between food sources.

## Usage

``` r
getEncounterByPrey(
  params,
  n = initialN(params),
  n_pp = initialNResource(params),
  n_other = initialNOther(params),
  t = 0,
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object.

- n:

  A matrix of species abundances (species x size). Defaults to the
  initial abundances stored in `params`.

- n_pp:

  A vector of the resource abundance by size. Defaults to the initial
  resource abundance stored in `params`.

- n_other:

  A named list of the abundances of other dynamical components. Defaults
  to the initial values stored in `params`.

- t:

  The time at which to evaluate the rates. Defaults to 0. It matters
  only in a model whose rates depend on time.

- ...:

  Passed on to the encounter functions of the components.

## Value

An array (predator species x predator size x food source) with the
encounter rates in grams per year. Dimnames are "predator", "w" and
"prey". The food sources are the prey species, "Resource", "External"
and the components that have an encounter rate function, in that order.
Unlike
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md), the
values are not set to zero outside the size range of the predator
species.

## Details

The rate \\E\_{ij}(w)\\ at which a predator of species \\i\\ and size
\\w\\ encounters biomass of prey species \\j\\ is \$\$ E\_{ij}(w) =
\gamma_i(w) \theta\_{ij} \int N_j(w_p) \phi_i(w, w_p) w_p dw_p, \$\$
with the notation of
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md),
and with \\\theta\_{ij}\\ replaced by the `interaction_resource` species
parameter for the resource. Further food sources are the external
encounter rate, set with
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
and each extra ecosystem component for which you have defined an
encounter rate function with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md).

## Extension hook

`getEncounterByPrey()` is an S3 generic. Its `MizerParams` method
resolves the standard encounter rate of
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md).
An extension package that changes the encounter rate with a
[`projectEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
method should give this generic a method that makes the same change and
calls [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html), so that
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md) and
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md)
can split the changed encounter rate between the food sources correctly.
Without such a method they find that the resolved rates no longer add up
to the encounter rate, scale all of them in proportion, and say so.
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md)
takes from this method only how much food each prey species and the
resource provide. Which prey sizes that food comes from it takes from
the predation kernel stored in the model.

This function does not itself check that its result adds up to
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md):
the `MizerParams` method is reached through
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) from extension
methods that have already changed the parameters, and the encounter rate
calculated from those would include the change twice.

## See also

[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)

Other summary functions:
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getMortByCause()`](https://sizespectrum.org/mizer/reference/getMortByCause.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
by_prey <- getEncounterByPrey(NS_params)
# Summed over food sources this is the encounter rate
all.equal(rowSums(by_prey, dims = 2), unclass(getEncounter(NS_params)),
          check.attributes = FALSE)
#> [1] TRUE
```
