# Get energy rate available for growth needed to project standard mizer model

Calculates the energy rate \\g_i(w)\\ (grams/year) available by species
and size for growth after metabolism, movement and reproduction have
been accounted for. Used by
[`project()`](https://sizespectrum.org/mizer/reference/project.md) for
performing simulations. You would not usually call this function
directly but instead use
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
which then calls this function unless an alternative function has been
registered, see below.

## Usage

``` r
projectEGrowth(params, n, n_pp, n_other, t = 0, e_repro, e, ...)

# S3 method for class 'MizerParams'
projectEGrowth(params, n, n_pp, n_other, t = 0, e_repro, e, ...)

mizerEGrowth(params, n, n_pp, n_other, t = 0, e_repro, e, ...)
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

- e_repro:

  The energy available for reproduction as calculated by
  [`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md).

- e:

  The energy available for reproduction and growth as calculated by
  [`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md).

- ...:

  Unused

## Value

A two dimensional array (species x size) with the growth rates.

## Details

The growth rate is calculated as the difference between the energy
available for reproduction and growth (obtainable with
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md))
and the energy used for reproduction (obtainable with
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md)),
but is set to 0 if the result would be negative.

## Your own growth rate function

By default
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md)
calls `mizerEGrowth()`. However you can replace this with your own
alternative growth rate function. If your function is called
`"myEGrowth"` then you register it in a MizerParams object `params` with

    params <- setRateFunction(params, "EGrowth", "myEGrowth")

Your function will then be called instead of `mizerEGrowth()`, with the
same arguments.

## See also

Other mizer rate functions:
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
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md),
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
