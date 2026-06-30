# Get energy rate available for reproduction and growth

Calculates the energy rate \\E\_{r.i}(w)\\ (grams/year) available for
reproduction and growth after metabolism and movement have been
accounted for.

## Usage

``` r
getEReproAndGrowth(object, ...)
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

## Value

- `MizerParams`: An `ArraySpeciesBySize` object (species x size) with
  the energy rate \\E\_{r.i}(w)\\ available for growth and reproduction
  (grams/year).

- `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species
  x size) with the energy rate at every time step. If `drop = TRUE` then
  dimensions of length 1 will be removed.

## Your own energy rate function

By default `getEReproAndGrowth()` calls
[`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md).
However you can replace this with your own alternative energy rate
function. If your function is called `"myEReproAndGrowth"` then you
register it in a MizerParams object `params` with

    params <- setRateFunction(params, "EReproAndGrowth", "myEReproAndGrowth")

Your function will then be called instead of
[`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md),
with the same arguments.

## See also

The part of this energy rate that is invested into growth is calculated
with
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md)
and the part that is invested into reproduction is calculated with
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md).

Other rate functions:
[`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md),
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
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
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
# \donttest{
params <- NS_params
# Project with constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the energy at a particular time step
e <- getEReproAndGrowth(params, n = N(sim)[15, , ],
                        n_pp = NResource(sim)[15, ], t = 15)
# Rate at this time for Sprat of size 2g
e["Sprat", "2"]
#> [1] 4.706336
# }
```
