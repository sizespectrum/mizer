# Get the total fishing mortality rate from all fishing gears

Calculates the total fishing mortality (in units 1/year) from all gears
by species and size. The total fishing mortality is just the sum of the
fishing mortalities imposed by each gear, \\\mu\_{f.i}(w)=\sum_g
F\_{g,i,w}\\. You would not usually call this function directly but
instead use
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md),
which then calls this function unless an alternative function has been
registered, see below.

## Usage

``` r
projectFMort(params, n, n_pp, n_other, t = 0, effort, e_growth, pred_mort, ...)

# S3 method for class 'MizerParams'
projectFMort(params, n, n_pp, n_other, t = 0, effort, e_growth, pred_mort, ...)

mizerFMort(params, n, n_pp, n_other, t = 0, effort, e_growth, pred_mort, ...)
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

- effort:

  A vector with the effort for each fishing gear.

- e_growth:

  An array (species x size) with the energy available for growth as
  calculated by
  [`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md).
  Unused.

- pred_mort:

  A two dimensional array (species x size) with the predation mortality
  as calculated by
  [`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md).
  Unused.

- ...:

  Unused

## Value

An array (species x size) with the fishing mortality.

## Note

Here: fishing mortality = catchability x selectivity x effort.

## Your own fishing mortality function

By default
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md)
calls `mizerFMort()`. However you can replace this with your own
alternative fishing mortality function. If your function is called
`"myFMort"` then you register it in a MizerParams object `params` with

    params <- setRateFunction(params, "FMort", "myFMort")

Your function will then be called instead of `mizerFMort()`, with the
same arguments.

## See also

Other mizer rate functions:
[`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md),
[`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md),
[`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md),
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md),
[`mizerFMortGear()`](https://sizespectrum.org/mizer/reference/mizerFMortGear.md),
[`mizerFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md),
[`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md),
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md),
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md),
[`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md),
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md),
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
