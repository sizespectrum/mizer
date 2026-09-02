# Superseded `get`-prefixed aliases for values stored in a model

**\[superseded\]**

Each of these functions is an alias for a function with a shorter name
that returns exactly the same thing. The shorter name is the one that
also has a replacement function (`catchability(params) <- value` and
friends), so that is the name to use:

|  |  |
|----|----|
| Superseded | Use instead |
| `getCatchability()` | [`catchability()`](https://sizespectrum.org/mizer/reference/setFishing.md) |
| `getSelectivity()` | [`selectivity()`](https://sizespectrum.org/mizer/reference/setFishing.md) |
| `getInitialEffort()` | [`initial_effort()`](https://sizespectrum.org/mizer/reference/initial_effort.md) |
| `getInteraction()` | [`interaction_matrix()`](https://sizespectrum.org/mizer/reference/setInteraction.md) |
| `getResourceDynamics()` | [`resource_dynamics()`](https://sizespectrum.org/mizer/reference/setResource.md) |
| `getResourceLevel()` | [`resource_level()`](https://sizespectrum.org/mizer/reference/setResource.md) |
| `getResourceRate()` | [`resource_rate()`](https://sizespectrum.org/mizer/reference/setResource.md) |
| `getResourceCapacity()` | [`resource_capacity()`](https://sizespectrum.org/mizer/reference/setResource.md) |
| `getPredKernel()` | [`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md) |
| `getSearchVolume()` | [`search_vol()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md) |
| `getMaxIntakeRate()` | [`intake_max()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md) |
| `getMetabolicRate()` | [`metab()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md) |
| `getExtMort()` | [`ext_mort()`](https://sizespectrum.org/mizer/reference/setExtMort.md) |
| `getExtEncounter()` | [`ext_encounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md) |
| `getMaturityProportion()` | [`maturity()`](https://sizespectrum.org/mizer/reference/setReproduction.md) |
| `getReproductionProportion()` | [`repro_prop()`](https://sizespectrum.org/mizer/reference/setReproduction.md) |
| `getReproductionLevel()` | [`reproduction_level()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md) |

The `get` prefix is reserved for the functions that *calculate* a rate
from the current state of a model, like
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
or [`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md).
The functions above only read back a value that is already stored in the
MizerParams object, which is what the bare names say.

The old names are however not going away. They are plain aliases: they
do not warn and they will keep working, so existing code and old scripts
run unchanged. They are not used anywhere inside mizer and should not be
used in new code.

## Usage

``` r
getCatchability(params)

getSelectivity(params)

getInitialEffort(params)

getInteraction(params)

getResourceDynamics(params)

getResourceLevel(params)

getResourceRate(params)

getResourceCapacity(params)

getPredKernel(params)

getSearchVolume(params)

getMaxIntakeRate(params)

getMetabolicRate(params)

getExtMort(params)

getExtEncounter(params)

getMaturityProportion(params)

getReproductionProportion(params)

getReproductionLevel(params)
```

## Arguments

- params:

  A MizerParams object

## Value

The same as the function it is an alias for.
