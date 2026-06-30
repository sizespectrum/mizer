# Calculate the SSB of species

Calculates the spawning stock biomass (SSB) for each species. For a
`MizerSim` object this is returned for every saved time; for a
`MizerParams` object it is calculated from the initial state. SSB is the
total mass of all mature individuals.

## Usage

``` r
getSSB(object)
```

## Arguments

- object:

  An object of class `MizerParams` or `MizerSim`.

## Value

If called with a MizerParams object, a named vector with the SSB in
grams for each species in the model. If called with a MizerSim object, a
`ArrayTimeBySpecies` object (time x species) containing the SSB in grams
at each time step for all species.

## See also

Other summary functions:
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
ssb <- getSSB(NS_sim)
ssb[c("1972", "2010"), c("Herring", "Cod")]
#> Spawning stock biomass (2 times x 2 species) [g] 
#>   Herring: min=4.43e+10 mean=1.09e+11 max=1.73e+11
#>   Cod: min=3.47e+11 mean=3.62e+11 max=3.76e+11
```
