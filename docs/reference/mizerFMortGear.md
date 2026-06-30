# Get the fishing mortality needed to project standard mizer model

Calculates the fishing mortality rate \\F\_{g,i,w}\\ by gear, species
and size. This is a helper function for
[`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md).

## Usage

``` r
mizerFMortGear(params, effort)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- effort:

  A vector with the effort for each fishing gear.

## Value

A three dimensional array (gear x species x size) with the fishing
mortality.

## Note

Here: fishing mortality = catchability x selectivity x effort.

## See also

[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)

Other mizer rate functions:
[`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md),
[`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md),
[`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md),
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md),
[`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md),
[`mizerFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md),
[`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md),
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md),
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md),
[`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md),
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md),
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
