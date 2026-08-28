# Description of summary functions

Mizer provides a range of functions to summarise the results of a
simulation.

## Details

A list of available summary functions is given in the table below.

|                                                                                                      |                                                                  |                                                                                                                                         |
|------------------------------------------------------------------------------------------------------|------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| Function                                                                                             | Returns                                                          | Description                                                                                                                             |
| [`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md)                                   | Three dimensional array (predator x size x prey)                 | Diet of predator at size, resolved by prey species                                                                                      |
| [`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md)                   | `ArraySpeciesBySize` (species x size)                            | Trophic level of individuals at size, accounting for ontogenetic diet shifts                                                            |
| [`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md) | Named vector (species)                                           | Consumption-rate-weighted mean trophic level of each species                                                                            |
| [`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md)                                     | Two dimensional array (time x species)                           | Total Spawning Stock Biomass (SSB) of each species through time where SSB is calculated as the sum of weight of all mature individuals. |
| [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)                             | Two dimensional array (time x species)                           | Total biomass of each species through time.                                                                                             |
| [`getN()`](https://sizespectrum.org/mizer/reference/getN.md)                                         | Two dimensional array (time x species)                           | Total abundance of each species through time.                                                                                           |
| [`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md)                   | Three dimensional array (time x species x size)                  | Feeding level of each species by size through time.                                                                                     |
| [`getM2`](https://sizespectrum.org/mizer/reference/getM2.md)                                         | Three dimensional array (time x species x size)                  | The predation mortality imposed on each species by size through time.                                                                   |
| [`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md)                                 | Three dimensional array (time x species x size)                  | Total fishing mortality on each species by size through time.                                                                           |
| [`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md)                         | Four dimensional array (time x gear x species x size)            | Fishing mortality on each species by each gear at size through time.                                                                    |
| [`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)                         | Three dimensional array (time x gear x species)                  | Total yield by gear and species through time.                                                                                           |
| [`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md)                                 | Two dimensional array (time x species)                           | Total yield of each species across all gears through time.                                                                              |
| [`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md)                         | Named vector (species) or two dimensional array (time x species) | Any integral over the size spectrum, from which all of the above are built. Use it to write your own summary function.                  |

## Writing your own summary function

The entry point is
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md).
It selects the size range, applies the bin-averaging appropriate to the
model's
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
setting and wraps the result in the right mizer array class, so a
summary function built on it is automatically consistent with the
quadrature the model is actually using. Pass the whole weighting factor
\\K(w)\\ evaluated on the size grid, but neither the bin widths
`params@dw` nor any bin-averaging of your own:
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md)
supplies both.

If your quantity involves the predation kernel, take the kernel from
[`encounter_kernel()`](https://sizespectrum.org/mizer/reference/encounter_kernel.md)
rather than from
[`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
and pair it with the plain point prey weight
`params@w_full * params@dw_full`. That weight is a normalisation which
the kernel construction is built to cancel, not a quadrature weight, so
it is the one place where you must *not* bin-average. Pairing the
point-sampled
[`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
with a bin-averaged prey weight applies the prey-bin integral twice.

## See also

[indicator_functions](https://sizespectrum.org/mizer/reference/indicator_functions.md),
[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.md)
