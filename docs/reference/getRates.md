# Get all rates

Calls other rate functions in sequence and collects the results in a
list. The rates returned are encounter, feeding level, energy for growth
and reproduction, predation rate, predation mortality, and resource
mortality. The purpose of this function is to provide a convenient way
to get all the rates at once, and to ensure that they are all calculated
at the same time step with the same inputs. The rates are returned in a
list with the same names as the rate functions that calculate them, so
for example the encounter rate is returned in the list element named
"encounter" and is calculated with the getEncounter() function.

## Usage

``` r
getRates(
  params,
  n = initialN(params),
  n_pp = initialNResource(params),
  n_other = initialNOther(params),
  effort,
  t = 0,
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components of the ecosystem

- effort:

  The effort for each fishing gear

- t:

  The time for which to do the calculation (Not used by standard mizer
  rate functions but useful for extensions with time-dependent
  parameters.)

- ...:

  Unused

## Details

When mizer needs to calculate the rates during a simulation it does not
use this function but instead the faster
[`projectRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md).

## See also

Other rate functions:
[`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md),
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md),
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md),
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md),
[`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md),
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md),
[`getFlux()`](https://sizespectrum.org/mizer/reference/getFlux.md),
[`getFluxGradient()`](https://sizespectrum.org/mizer/reference/getFluxGradient.md),
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md),
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md),
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md),
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md),
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
rates <- getRates(NS_params)
names(rates)
#>  [1] "encounter"     "feeding_level" "e"             "e_repro"      
#>  [5] "e_growth"      "diffusion"     "pred_rate"     "pred_mort"    
#>  [9] "f_mort"        "mort"          "rdi"           "rdd"          
#> [13] "resource_mort"
identical(rates$encounter, getEncounter(NS_params))
#> [1] FALSE
```
