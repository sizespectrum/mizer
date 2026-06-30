# Alias for `getResourceMort()`

**\[deprecated\]** An alias provided for backward compatibility with
mizer version \<= 1.0

## Usage

``` r
getM2Background(
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
  object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components of the ecosystem

- t:

  The time for which to do the calculation (Not used by standard mizer
  rate functions but useful for extensions with time-dependent
  parameters.)

- ...:

  Unused

## Value

A vector of mortality rate by resource size.

## Your own resource mortality function

By default
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)
calls
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md).
However you can replace this with your own alternative resource
mortality function. If your function is called `"myResourceMort"` then
you register it in a MizerParams object `params` with

    params <- setRateFunction(params, "ResourceMort", "myResourceMort")

Your function will then be called instead of
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md),
with the same arguments.

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
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md)

## Examples

``` r
# \donttest{
params <- NS_params
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get resource mortality at one time step
getResourceMort(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ])
#> Resource mortality (226 sizes) [1/year] 
#>   min=3.9e-16 mean=5.5 max=17.6
# }
```
