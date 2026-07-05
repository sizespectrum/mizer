# Calculate the rate at which biomass of each species is fished by each gear

This yield rate is given in grams per year. It is calculated at each
time step saved in the MizerSim object.

## Usage

``` r
getYieldGear(object)
```

## Arguments

- object:

  An object of class `MizerParams` or `MizerSim`.

## Value

If called with a MizerParams object, an array (gear x species) with the
yield rate in grams per year from each gear for each species in the
model. If called with a MizerSim object, an array (time x gear x
species) containing the yield rate at each time step.

## Details

For details of how the yield rate is defined see the help page of
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md).

## See also

[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md)

Other summary functions:
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md)

## Examples

``` r
yield <- getYieldGear(NS_sim)
dim(yield)
#> [1] 44 12 12
yield["1972", , "Herring"]
#>       Sprat     Sandeel      N.pout     Herring         Dab     Whiting 
#>           0           0           0 80002050975           0           0 
#>        Sole     Gurnard      Plaice     Haddock         Cod      Saithe 
#>           0           0           0           0           0           0 
```
