# Get total predation mortality rate

Calculates the total predation mortality rate \\\mu\_{p,i}(w_p)\\ (in
units of 1/year) on each prey species by prey size: \$\$\mu\_{p.i}(w_p)
= \sum_j {\tt pred\\rate}\_j(w_p)\\ \theta\_{ji}.\$\$ The predation rate
`pred_rate` is returned by
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md).

## Usage

``` r
getPredMort(object, ...)
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

- `MizerParams`: An `ArraySpeciesBySize` object (prey species x prey
  size) with the predation mortality rates.

- `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x prey
  species x prey size) with the predation mortality at every time step.
  If `drop = TRUE` then dimensions of length 1 will be removed.

## Your own predation mortality function

By default `getPredMort()` calls
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md).
However you can replace this with your own alternative predation
mortality function. If your function is called `"myPredMort"` then you
register it in a MizerParams object `params` with

    params <- setRateFunction(params, "PredMort", "myPredMort")

Your function will then be called instead of
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md),
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
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md),
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md),
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md),
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
# \donttest{
params <- NS_params
# Predation mortality in initial state
M2 <- getPredMort(params)
str(M2)
#>  'ArraySpeciesBySize' num [1:12, 1:100] 3.64 4.43 4.31 4.89 4.8 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ sp: chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>  - attr(*, "value_name")= chr "Predation mortality"
#>  - attr(*, "units")= chr "1/year"
#>  - attr(*, "representation")= chr "average"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get predation mortality at one time step
M2 <- getPredMort(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ])
# Get predation mortality at all saved time steps
M2 <- getPredMort(sim)
str(M2)
#>  'ArrayTimeBySpeciesBySize' num [1:21, 1:12, 1:100] 3.64 3.42 3.29 3.2 3.16 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ time: chr [1:21] "0" "1" "2" "3" ...
#>   ..$ sp  : chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w   : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>  - attr(*, "value_name")= chr "Predation mortality"
#>  - attr(*, "units")= chr "1/year"
#>  - attr(*, "representation")= chr "average"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
# Get predation mortality over the years 15 - 20
M2 <- getPredMort(sim, time_range = c(15, 20))
# }
```
