# Get the total fishing mortality rate from all fishing gears by time, species and size.

Calculates the total fishing mortality (in units 1/year) from all gears
by species and size and possibly time. See
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
for details of how fishing gears are set up.

## Usage

``` r
getFMort(object, ...)
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

  :   The effort of each fishing gear. See notes below. Defaults to the
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

  `drop`

  :   Should dimensions of length 1 be dropped, e.g. if your community
      only has one species it might make presentation of results easier.
      Defaults to `TRUE`.

## Value

- `MizerParams` with vector effort: An `ArraySpeciesBySize` object
  (species x size) with the fishing mortality rates.

- `MizerParams` with time-dimensioned effort or `MizerSim`: An
  `ArrayTimeBySpeciesBySize` object (time x species x size).

The `effort` argument is only used if a `MizerParams` object is passed
in. The `effort` argument can be a two dimensional array (time x gear),
a vector of length equal to the number of gears (each gear has a
different effort that is constant in time), or a single numeric value
(each gear has the same effort that is constant in time). The order of
gears in the `effort` argument must be the same as in the `MizerParams`
object.

If the object argument is of class `MizerSim` then the effort slot of
the `MizerSim` object is used and the `effort` argument is not used.

## Details

The total fishing mortality is just the sum of the fishing mortalities
imposed by each gear, \\F_i(w)=\sum_g F\_{g,i,w}\\. The fishing
mortality for each gear is obtained as catchability x selectivity x
effort.

## Your own fishing mortality function

By default `getFMort()` calls
[`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md).
However you can replace this with your own alternative fishing mortality
function. If your function is called `"myFMort"` then you register it in
a MizerParams object `params` with

    params <- setRateFunction(params, "FMort", "myFMort")

Your function will then be called instead of
[`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md),
with the same arguments.

## See also

Other rate functions:
[`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md),
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md),
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md),
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
# Get the total fishing mortality in the initial state
F <- getFMort(params, effort = 1)
str(F)
#>  'ArraySpeciesBySize' num [1:12, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ sp: chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>  - attr(*, "value_name")= chr "Fishing mortality"
#>  - attr(*, "units")= chr "1/year"
#>  - attr(*, "representation")= chr "average"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
# Get the initial total fishing mortality when effort is different
# between the four gears:
F <- getFMort(params, effort = c(0.5,1,1.5,0.75))
# Get the total fishing mortality when effort is different
# between the four gears and changes with time:
effort <- array(NA, dim = c(20,4))
effort[, 1] <- seq(from = 0, to = 1, length = 20)
effort[, 2] <- seq(from = 1, to = 0.5, length = 20)
effort[, 3] <- seq(from = 1, to = 2, length = 20)
effort[, 4] <- seq(from = 2, to = 1, length = 20)
F <- getFMort(params, effort = effort)
str(F)
#>  num [1:20, 1:12, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ time: NULL
#>   ..$ sp  : chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w   : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
# Get the total fishing mortality using the effort already held in a
# MizerSim object.
sim <- project(params, t_max = 20, effort = 0.5)
F <- getFMort(sim)
F <- getFMort(sim, time_range = c(10, 20))
# }
```
