# Get density independent rate of egg production

Calculates the density-independent rate of total egg production
\\R\_{di}\\ (units 1/year) before density dependence, by species.

## Usage

``` r
getRDI(object, ...)
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

## Value

- `MizerParams`: A numeric vector the length of the number of species.

- `MizerSim`: An `ArrayTimeBySpecies` object (time x species).

## Details

This rate is obtained by taking the per capita rate \\E_r(w)\psi(w)\\ at
which energy is invested in reproduction, as calculated by
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
multiplying it by the number of individuals\\N(w)\\ and integrating over
all sizes \\w\\ and then multiplying by the reproductive efficiency
\\\epsilon\\ and dividing by the egg size `w_min`, and by a factor of
two to account for the two sexes: \$\$R\_{di} = \frac{\epsilon}{2
w\_{min}} \int N(w) E_r(w) \psi(w) \\ dw\$\$

Used by [`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md)
to calculate the actual, density dependent rate. See
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
for more details.

## Your own reproduction function

By default `getRDI()` calls
[`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md).
However you can replace this with your own alternative reproduction
function. If your function is called `"myRDI"` then you register it in a
MizerParams object `params` with

    params <- setRateFunction(params, "RDI", "myRDI")

Your function will then be called instead of
[`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md),
with the same arguments. For an example of an alternative reproduction
function see
[`constantEggRDI()`](https://sizespectrum.org/mizer/reference/constantEggRDI.md).

## See also

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
[`getFlux()`](https://sizespectrum.org/mizer/reference/getFlux.md),
[`getFluxGradient()`](https://sizespectrum.org/mizer/reference/getFluxGradient.md),
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md),
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md),
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md),
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md),
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
# \donttest{
params <- NS_params
# Project with constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the density-independent reproduction rate at a particular time step
getRDI(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
#>        Sprat      Sandeel       N.pout      Herring          Dab      Whiting 
#> 5.041129e+13 1.053968e+15 9.020326e+13 2.024648e+14 2.517323e+12 3.476507e+13 
#>         Sole      Gurnard       Plaice      Haddock          Cod       Saithe 
#> 9.988220e+12 1.387214e+12 2.839019e+13 2.983237e+13 1.143392e+14 4.340615e+13 
# }
```
