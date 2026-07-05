# Description of summary functions

Mizer provides a range of functions to summarise the results of a
simulation.

## Details

A list of available summary functions is given in the table below.

|  |  |  |
|----|----|----|
| Function | Returns | Description |
| [`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md) | Three dimensional array (predator x size x prey) | Diet of predator at size, resolved by prey species |
| [`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md) | `ArraySpeciesBySize` (species x size) | Trophic level of individuals at size, accounting for ontogenetic diet shifts |
| [`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md) | Named vector (species) | Consumption-rate-weighted mean trophic level of each species |
| [`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md) | Two dimensional array (time x species) | Total Spawning Stock Biomass (SSB) of each species through time where SSB is calculated as the sum of weight of all mature individuals. |
| [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md) | Two dimensional array (time x species) | Total biomass of each species through time. |
| [`getN()`](https://sizespectrum.org/mizer/reference/getN.md) | Two dimensional array (time x species) | Total abundance of each species through time. |
| [`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md) | Three dimensional array (time x species x size) | Feeding level of each species by size through time. |
| [`getM2`](https://sizespectrum.org/mizer/reference/getM2.md) | Three dimensional array (time x species x size) | The predation mortality imposed on each species by size through time. |
| [`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md) | Three dimensional array (time x species x size) | Total fishing mortality on each species by size through time. |
| [`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md) | Four dimensional array (time x gear x species x size) | Fishing mortality on each species by each gear at size through time. |
| [`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md) | Three dimensional array (time x gear x species) | Total yield by gear and species through time. |
| [`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md) | Two dimensional array (time x species) | Total yield of each species across all gears through time. |

## See also

[indicator_functions](https://sizespectrum.org/mizer/reference/indicator_functions.md),
[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.md)
