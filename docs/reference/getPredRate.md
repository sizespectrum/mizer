# Get predation rate

Calculates the potential rate (in units 1/year) at which a prey
individual of a given size \\w\\ is killed by predators from species
\\j\\. In formulas \$\${\tt pred\\rate}\_j(w_p) = \int \phi_j(w,w_p)
(1-f_j(w)) \gamma_j(w) N_j(w) \\ dw.\$\$ This potential rate is used in
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md)
to calculate the realised predation mortality rate on the prey
individual.

## Usage

``` r
getPredRate(object, ...)
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

- `MizerParams`: An `ArraySpeciesBySize` object (predator species x prey
  size), where the prey size runs over fish community plus resource
  spectrum.

- `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x predator
  species x prey size) with the predation rates at every time step. If
  `drop = TRUE` then dimensions of length 1 will be removed.

## Your own predation rate function

By default `getPredRate()` calls
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md).
However you can replace this with your own alternative predation rate
function. If your function is called `"myPredRate"` then you register it
in a MizerParams object `params` with

    params <- setRateFunction(params, "PredRate", "myPredRate")

Your function will then be called instead of
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md),
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
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md),
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md),
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
# \donttest{
params <- NS_params
# Predation rate in initial state
pred_rate <- getPredRate(params)
str(pred_rate)
#>  'ArraySpeciesBySize' num [1:12, 1:226] 8.35e-17 6.05e-10 9.75e-16 1.19e-05 1.04e-17 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ sp    : chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w_prey: chr [1:226] "2.12e-13" "2.53e-13" "3.02e-13" "3.61e-13" ...
#>  - attr(*, "value_name")= chr "Predation rate"
#>  - attr(*, "units")= chr "1/year"
#>  - attr(*, "representation")= chr "point"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the feeding level at one time step
pred_rate <- getPredRate(params, n = N(sim)[15, , ],
                         n_pp = NResource(sim)[15, ], t = 15)
# }
```
