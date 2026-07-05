# Get predation mortality rate for resource needed to project standard mizer model

Calculates the predation mortality rate \\\mu_p(w)\\ on the resource
spectrum by resource size (in units 1/year). You would not usually call
this function directly but instead use
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md),
which then calls this function unless an alternative function has been
registered, see below.

## Usage

``` r
projectResourceMort(params, n, n_pp, n_other, t = 0, pred_rate, ...)

# S3 method for class 'MizerParams'
projectResourceMort(params, n, n_pp, n_other, t = 0, pred_rate, ...)

mizerResourceMort(params, n, n_pp, n_other, t = 0, pred_rate, ...)
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

- pred_rate:

  A two dimensional array (predator species x prey size) with the
  predation rate, where the prey size runs over fish community plus
  resource spectrum.

- ...:

  Unused

## Value

A vector of mortality rate by resource size.

## Your own resource mortality function

By default
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)
calls `mizerResourceMort()`. However you can replace this with your own
alternative resource mortality function. If your function is called
`"myResourceMort"` then you register it in a MizerParams object `params`
with

    params <- setRateFunction(params, "ResourceMort", "myResourceMort")

Your function will then be called instead of `mizerResourceMort()`, with
the same arguments.

## See also

Other mizer rate functions:
[`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md),
[`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md),
[`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md),
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md),
[`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md),
[`mizerFMortGear()`](https://sizespectrum.org/mizer/reference/mizerFMortGear.md),
[`mizerFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md),
[`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md),
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md),
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md),
[`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md),
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md)
