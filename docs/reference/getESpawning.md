# Alias for `getERepro()`

**\[deprecated\]** An alias provided for backward compatibility with
mizer version \<= 1.0

## Usage

``` r
getESpawning(object, ...)
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

- `MizerParams`: An `ArraySpeciesBySize` object (species x size) holding
  \$\$\psi_i(w)\max(0, E\_{r.i}(w))\$\$ where \\E\_{r.i}(w)\\ is the
  rate at which energy becomes available for growth and reproduction,
  calculated with
  [`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md),
  and \\\psi_i(w)\\ is the proportion of this energy that is used for
  reproduction. Negative values of \\E\_{r.i}(w)\\ are clipped to 0
  before multiplying by \\\psi_i(w)\\. This proportion is taken from the
  `params` object and is set with
  [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md).

- `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species
  x size) with the energy for reproduction at every time step. If
  `drop = TRUE` then dimensions of length 1 will be removed.

## Your own reproduction rate function

By default
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md)
calls
[`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md).
However you can replace this with your own alternative reproduction rate
function. If your function is called `"myERepro"` then you register it
in a MizerParams object `params` with

    params <- setRateFunction(params, "ERepro", "myERepro")

Your function will then be called instead of
[`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md),
with the same arguments.

## See also

Other rate functions:
[`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md),
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
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
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
# \donttest{
params <- NS_params
# Project with constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the rate at a particular time step
erepro <- getERepro(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
# Rate at this time for Sprat of size 2g
erepro["Sprat", "2"]
#> [1] 0
# }
```
