# Get total mortality rate

Calculates the total mortality rate \\\mu_i(w)\\ (in units 1/year) on
each species by size from predation mortality, background mortality and
fishing mortality for a single time step.

## Usage

``` r
getMort(object, ...)
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

  `effort`

  :   A numeric vector of the effort by gear or a single numeric effort
      value which is used for all gears. Defaults to the initial effort
      stored in `object`.

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
  the total mortality rates.

- `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species
  x size) with the total mortality rates at every time step. If
  `drop = TRUE` then dimensions of length 1 will be removed.

## Details

If your model contains additional components that you added with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
and for which you specified a `mort_fun` function then the mortality
inflicted by these components will be included in the returned value.

## Your own mortality function

By default `getMort()` calls
[`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md).
However you can replace this with your own alternative mortality
function. If your function is called `"myMort"` then you register it in
a MizerParams object `params` with

    params <- setRateFunction(params, "Mort", "myMort")

Your function will then be called instead of
[`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md),
with the same arguments.

## See also

[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md),
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md)

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
# Get the total mortality at a particular time step
mort <- getMort(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ],
                t = 15, effort = 0.5)
# Mortality rate at this time for Sprat of size 2g
mort["Sprat", "2"]
#> [1] 1.640475
# }
```
