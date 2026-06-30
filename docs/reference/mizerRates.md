# Get all rates needed to project standard mizer model

Calls other rate functions in sequence and collects the results in a
list.

`projectRates()` is an S3 generic used by extension-aware projections to
calculate all rates. Models without extensions keep using `mizerRates()`
directly. The base method mirrors `mizerRates()` but calls migrated
projection hooks directly, starting with
[`projectEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md).

## Usage

``` r
mizerRates(params, n, n_pp, n_other, t = 0, effort, rates_fns, ...)

projectRates(params, n, n_pp, n_other, t = 0, effort, rates_fns, ...)
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

  The effort for each fishing gear

- rates_fns:

  Named list of the functions to call to calculate the rates. Note that
  this list holds the functions themselves, not their names.

- ...:

  Unused

## Value

List of rates.

## Details

By default this function returns a list with the following components:

- encounter from
  [`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)

- feeding_level from
  [`mizerFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md)

- e from
  [`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md)

- e_repro from
  [`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md)

- e_growth from
  [`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md)

- pred_rate from
  [`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md)

- pred_mort from
  [`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md)

- f_mort from
  [`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md)

- mort from
  [`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md)

- rdi from
  [`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md)

- rdd from
  [`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md)

- resource_mort from
  [`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)

However you can replace any of these rate functions by your own rate
function if you wish, see
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
for details.

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
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
