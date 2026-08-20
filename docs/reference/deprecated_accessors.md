# Deprecated `get`-prefixed accessors for stored rate arrays

**\[deprecated\]**

Each of these functions had a shorter alias returning exactly the same
thing. The shorter name is the one that also has a replacement function
(`catchability(params) <- value` and friends), so that is the name that
survives. Use the replacement given below instead:

|  |  |
|----|----|
| Deprecated | Use instead |
| `getCatchability()` | [`catchability()`](https://sizespectrum.org/mizer/reference/setFishing.md) |
| `getSelectivity()` | [`selectivity()`](https://sizespectrum.org/mizer/reference/setFishing.md) |
| `getInitialEffort()` | [`initial_effort()`](https://sizespectrum.org/mizer/reference/initial_effort.md) |
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
The functions above only read back a value stored in the MizerParams
object, which is what the bare names say.

## Usage

``` r
getCatchability(params)

getSelectivity(params)

getInitialEffort(params)

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

The same as the replacement function it forwards to.
