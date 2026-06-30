# Get density dependent reproduction rate

Calculates the density dependent rate of egg production \\R_i\\ (units
1/year) for each species. This is the flux entering the smallest size
class of each species. The density dependent rate is the density
independent rate obtained with
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md) after
it has been put through the density dependence function. This is the
Beverton-Holt function
[`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md)
by default, but this can be changed. See
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
for more details.

## Usage

``` r
getRDD(object, ...)
```

## Arguments

- object:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  or
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  object.

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

  `t`

  :   The time for which to do the calculation. Defaults to 0.

  `rdi`

  :   A vector of density-independent reproduction rates for each
      species. If not specified, it is calculated internally using
      [`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md).

  **For a
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  object:**

  `time_range`

  :   The time range over which to return the rates. Either a vector of
      values, a vector of min and max time, or a single value. Defaults
      to the whole time range of the simulation.

## Value

- `MizerParams`: A numeric vector the length of the number of species.

- `MizerSim`: An `ArrayTimeBySpecies` object (time x species).

## See also

[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md)

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
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md),
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
# \donttest{
params <- NS_params
# Project with constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the rate at a particular time step
getRDD(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
#>        Sprat      Sandeel       N.pout      Herring          Dab      Whiting 
#> 7.273519e+11 4.098406e+11 9.405199e+12 1.103948e+12 1.115039e+10 5.394960e+11 
#>         Sole      Gurnard       Plaice      Haddock          Cod       Saithe 
#> 3.855063e+10 7.536194e+11 2.654321e+13 1.733106e+12 8.259403e+09 1.117118e+11 
# }
```
