# Get the fishing mortality by time, gear, species and size

Calculates the fishing mortality rate \\F\_{g,i,w}\\ by gear, species
and size and possibly time (in units 1/year).

## Usage

``` r
getFMortGear(object, ...)
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

  `effort`

  :   The effort for each fishing gear. See notes below. Defaults to the
      initial effort stored in `object`.

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

  :   Subset the returned fishing mortalities by time. The time range is
      either a vector of values, a vector of min and max time, or a
      single value. Defaults to the whole time range of the simulation.

## Value

An array. If the effort argument has a time dimension, or a `MizerSim`
is passed in, the output array has four dimensions (time x gear x
species x size). If the effort argument does not have a time dimension
(i.e. it is a vector or a single numeric), the output array has three
dimensions (gear x species x size).

## Note

Here: fishing mortality = catchability x selectivity x effort.

The `effort` argument is only used if a `MizerParams` object is passed
in. The `effort` argument can be a two dimensional array (time x gear),
a vector of length equal to the number of gears (each gear has a
different effort that is constant in time), or a single numeric value
(each gear has the same effort that is constant in time). The order of
gears in the `effort` argument must be the same the same as in the
`MizerParams` object. If the `effort` argument is not supplied, its
value is taken from the `@initial_effort` slot in the params object.

If the object argument is of class `MizerSim` then the effort slot of
the `MizerSim` object is used and the `effort` argument is not used.

## See also

Other rate functions:
[`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md),
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md),
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md),
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md),
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
params <-NS_params
# Get the fishing mortality in initial state
F <- getFMortGear(params, effort = 1)
str(F)
#>  num [1:4, 1:12, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ gear: chr [1:4] "Industrial" "Pelagic" "Beam" "Otter"
#>   ..$ sp  : chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w   : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
# Get the initial fishing mortality when effort is different
# between the four gears:
F <- getFMortGear(params, effort = c(0.5, 1, 1.5, 0.75))
# Get the fishing mortality when effort is different
# between the four gears and changes with time:
effort <- array(NA, dim = c(20, 4))
effort[, 1] <- seq(from=0, to = 1, length = 20)
effort[, 2] <- seq(from=1, to = 0.5, length = 20)
effort[, 3] <- seq(from=1, to = 2, length = 20)
effort[, 4] <- seq(from=2, to = 1, length = 20)
F <- getFMortGear(params, effort = effort)
str(F)
#>  num [1:20, 1:4, 1:12, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>  - attr(*, "dimnames")=List of 4
#>   ..$ time: NULL
#>   ..$ gear: chr [1:4] "Industrial" "Pelagic" "Beam" "Otter"
#>   ..$ sp  : chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w   : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
# Get the fishing mortality using the effort already held in a MizerSim object.
sim <- project(params, t_max = 20, effort = 0.5)
F <- getFMortGear(sim)
F <- getFMortGear(sim, time_range = c(10, 20))
# }
```
