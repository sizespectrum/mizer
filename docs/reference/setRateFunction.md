# Set own rate function to replace mizer rate function

If the way mizer calculates a fundamental rate entering the model is not
flexible enough for you (for example if you need to introduce time
dependence) then you can write your own functions for calculating that
rate and use `setRateFunction()` to register it with mizer.

## Usage

``` r
setRateFunction(params, rate, fun)

getRateFunction(params, rate)

other_params(params)

other_params(params) <- value
```

## Arguments

- params:

  A MizerParams object

- rate:

  Name of the rate for which a new function is to be set.

- fun:

  Name of the function to use to calculate the rate.

- value:

  A named list of user-defined parameters to store in
  `other_params(params)`.

## Value

For `setRateFunction()`: An updated MizerParams object

For `getRateFunction()`: The name of the registered rate function for
the requested `rate`, or the list of all rate functions if called
without `rate` argument.

For `other_params()`: The user-defined parameters stored in
`other_params(params)`, or `NULL` if none have been set. This excludes
any component-specific parameters stored via
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md).

## Details

At each time step during a simulation with the
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
function, mizer needs to calculate the instantaneous values of the
various rates. By default it calls the
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md)
function which creates a list with the following components:

- `encounter` from
  [`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)

- `feeding_level` from
  [`mizerFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md)

- `pred_rate` from
  [`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md)

- `pred_mort` from
  [`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md)

- `f_mort` from
  [`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md)

- `mort` from
  [`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md)

- `resource_mort` from
  [`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)

- `e` from
  [`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md)

- `e_repro` from
  [`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md)

- `e_growth` from
  [`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md)

- `diffusion` from
  [`mizerDiffusion()`](https://sizespectrum.org/mizer/reference/mizerDiffusion.md)

- `rdi` from
  [`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md)

- `rdd` from
  [`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md)

For each of these you can substitute your own function. So for example
if you have written your own function for calculating the total
mortality rate and have called it `myMort` and have a mizer model stored
in a MizerParams object called `params` that you want to run with your
new mortality rate, then you would call

    params <- setRateFunction(params, "Mort", "myMort")

In general if you want to replace a function `mizerSomeRateFunc()` with
a function `myVersionOfThis()` you would call

    params <- setRateFunction(params, "SomeRateFunc", "myVersionOfThis")

In some extreme cases you may need to swap out the entire
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md)
function for your own function called `myRates()`. That you can do with

    params <- setRateFunction(params, "Rates", "myRates")

Your new rate functions may need their own model parameters. These you
can store in `other_params(params)`. For example

    other_params(params)$my_param <- 42

Note that your own rate functions need to be defined in the global
environment or in a package. If they are defined within a function then
mizer will not find them.

## Avoid rates that jump as a function of abundance

Make sure your rate function depends **continuously** on the abundances
`n`, `n_pp` and `n_other`. It is tempting to write a rate that switches
abruptly on the state of the model — a fishery that closes when a stock
falls below a limit reference point, a predator that switches diet when
its preferred prey becomes scarce, a mortality that kicks in below a
critical condition. Such a rate breaks the assumption underlying every
one of mizer's time-stepping methods, which freeze the rates during each
density update and so cannot see a threshold being crossed within a
step.

The symptoms are quiet: mizer issues no warning, but the trajectory
keeps changing as you refine `dt`, the Newton solver stalls, and
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
reports a confident but meaningless answer. Choosing the L-stable
`method = "tr_bdf2"` in
[`project()`](https://sizespectrum.org/mizer/reference/project.md) does
not help, because the difficulty lies in the frozen rates rather than in
the linear solve.

The remedy is to give the switch a finite width, using a linear ramp
between two thresholds or a logistic transition, which is usually the
more realistic model anyway. Rates built with
[`max()`](https://rdrr.io/r/base/Extremes.html) or
[`min()`](https://rdrr.io/r/base/Extremes.html) are continuous but not
differentiable; these are much less troublesome, costing some accuracy
but not correctness. See the [Discontinuous rate
functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html)
article for the full story, the diagnostics and the fix.

## See also

"Extending mizer": [guide to extending
mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.html);
[Discontinuous rate
functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`clearExtensionChain()`](https://sizespectrum.org/mizer/reference/clearExtensionChain.md),
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`other_mort()`](https://sizespectrum.org/mizer/reference/other_mort.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md),
[`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
