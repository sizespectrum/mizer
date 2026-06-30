# Get flux into size bins

Calculates the flux \\J_i(w)\\ (numbers/year) entering each size class
from the one below it. This is composed of an advective flux from
somatic growth and a diffusive flux from the redistribution of
individuals.

## Usage

``` r
getFlux(object, ..., power = 0)
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

  **For a
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  object:**

  `time_range`

  :   The time range over which to return the rates. Either a vector of
      values, a vector of min and max time, or a single value. Defaults
      to the whole time range of the simulation.

  `drop`

  :   If `TRUE` then any dimension of length 1 is removed from the
      returned array.

- power:

  The flux at weight \\w\\ is multiplied by \\w\\ raised to `power`. The
  default `power = 0` gives the flux of individuals (numbers/year),
  whereas `power = 1` gives the flux of biomass (grams/year).

## Value

- `MizerParams`: An `ArraySpeciesBySize` object (species x size) with
  the flux entering each size class. The units are `numbers/year` when
  `power = 0` and `g^power/year` otherwise.

- `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species
  x size) with the flux at every time step. If `drop = TRUE` then
  dimensions of length 1 will be removed.

## Details

At the recruitment size, the flux is simply the recruitment rate
\\R\_{dd,i}\\ (see
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md)). For
sizes below the recruitment size the flux is zero.

The flux at weight \\w\\ is multiplied by \\w\\ raised to the `power`
given by the `power` argument, similar to the `power` argument of
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md).
The default `power = 0` returns the flux of individuals (numbers/year).
With `power = 1` the result is the flux of biomass (grams/year).

## See also

[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md)

Other rate functions:
[`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md),
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md),
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md),
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md),
[`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md),
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md),
[`getFluxGradient()`](https://sizespectrum.org/mizer/reference/getFluxGradient.md),
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md),
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md),
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md),
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md),
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md),
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
# \donttest{
params <- NS_params
# Project with constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the flux at a particular time step
flux <- getFlux(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
# Flux for Sprat of size 2g
flux["Sprat", "2"]
#> [1] 45461076783
# }
```
