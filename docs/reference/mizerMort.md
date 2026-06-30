# Get total mortality rate needed to project standard mizer model

Calculates the total mortality rate \\\mu_i(w)\\ (in units 1/year) on
each species by size from predation mortality, background mortality and
fishing mortality. You would not usually call this function directly but
instead use
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md),
which then calls this function unless an alternative function has been
registered, see below.

## Usage

``` r
projectMort(params, n, n_pp, n_other, t = 0, f_mort, pred_mort, ...)

# S3 method for class 'MizerParams'
projectMort(params, n, n_pp, n_other, t = 0, f_mort, pred_mort, ...)

mizerMort(params, n, n_pp, n_other, t = 0, f_mort, pred_mort, ...)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components of the ecosystem

- t:

  The time for which to do the calculation (Not used by standard mizer
  rate functions but useful for extensions with time-dependent
  parameters.)

- f_mort:

  A two dimensional array (species x size) with the fishing mortality

- pred_mort:

  A two dimensional array (species x size) with the predation mortality

- ...:

  Unused

## Value

A named two dimensional array (species x size) with the total mortality
rates.

## Details

If your model contains additional components that you added with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
and for which you specified a `mort_fun` function then the mortality
inflicted by these components will be included in the returned value.

## Your own mortality function

By default
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md) calls
`mizerMort()`. However you can replace this with your own alternative
mortality function. If your function is called `"myMort"` then you
register it in a MizerParams object `params` with

    params <- setRateFunction(params, "Mort", "myMort")

Your function will then be called instead of `mizerMort()`, with the
same arguments.

## See also

Other mizer rate functions:
[`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md),
[`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md),
[`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md),
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md),
[`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md),
[`mizerFMortGear()`](https://sizespectrum.org/mizer/reference/mizerFMortGear.md),
[`mizerFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md),
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md),
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md),
[`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md),
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md),
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
