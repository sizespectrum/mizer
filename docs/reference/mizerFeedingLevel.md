# Get feeding level needed to project standard mizer model

You would not usually call this function directly but instead use
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md),
which then calls this function unless an alternative function has been
registered, see below.

## Usage

``` r
projectFeedingLevel(params, n, n_pp, n_other, t = 0, encounter, ...)

# S3 method for class 'MizerParams'
projectFeedingLevel(params, n, n_pp, n_other, t = 0, encounter, ...)

mizerFeedingLevel(params, n, n_pp, n_other, t = 0, encounter, ...)
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

- encounter:

  A two dimensional array (predator species x predator size) with the
  encounter rate.

- ...:

  Unused

## Value

A two dimensional array (predator species x predator size) with the
feeding level.

## Feeding level

The feeding level \\f_i(w)\\ is the proportion of its maximum intake
rate at which the predator is actually taking in fish. It is calculated
from the encounter rate \\E_i\\ and the maximum intake rate \\h_i(w)\\
as \$\$f_i(w) = \frac{E_i(w)}{E_i(w)+h_i(w)}.\$\$ The encounter rate
\\E_i\\ is passed as an argument or calculated with
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md).
The maximum intake rate \\h_i(w)\\ is taken from the `params` object,
and is set with
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md).
As a consequence of the above expression for the feeding level,
\\1-f_i(w)\\ is the proportion of the food available to it that the
predator actually consumes.

## Your own feeding level function

By default
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md)
calls `mizerFeedingLevel()`. However you can replace this with your own
alternative feeding level function. If your function is called
`"myFeedingLevel"` then you register it in a MizerParams object `params`
with

    params <- setRateFunction(params, "FeedingLevel", "myFeedingLevel")

Your function will then be called instead of `mizerFeedingLevel()`, with
the same arguments.

## See also

The feeding level is used in
[`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md)
and in
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md).

Other mizer rate functions:
[`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md),
[`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md),
[`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md),
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md),
[`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md),
[`mizerFMortGear()`](https://sizespectrum.org/mizer/reference/mizerFMortGear.md),
[`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md),
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md),
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md),
[`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md),
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md),
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
