# Get energy rate available for reproduction and growth needed to project standard mizer model

Calculates the energy rate \\E\_{r.i}(w)\\ (grams/year) available to an
individual of species i and size w for reproduction and growth after
metabolism and movement have been accounted for. You would not usually
call this function directly but instead use
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md),
which then calls this function unless an alternative function has been
registered, see below.

## Usage

``` r
projectEReproAndGrowth(
  params,
  n,
  n_pp,
  n_other,
  t = 0,
  encounter,
  feeding_level,
  ...
)

# S3 method for class 'MizerParams'
projectEReproAndGrowth(
  params,
  n,
  n_pp,
  n_other,
  t = 0,
  encounter,
  feeding_level,
  ...
)

mizerEReproAndGrowth(
  params,
  n,
  n_pp,
  n_other,
  t = 0,
  encounter,
  feeding_level,
  ...
)
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

  An array (species x size) with the encounter rate as calculated by
  [`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md).

- feeding_level:

  An array (species x size) with the feeding level as calculated by
  [`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md).

- ...:

  Unused

## Value

A two dimensional array (species x size) holding \$\$E\_{r.i}(w) =
\alpha_i\\ (1 - {\tt feeding\\level}\_i(w))\\ {\tt encounter}\_i(w) -
{\tt metab}\_i(w).\$\$ Due to the form of the feeding level, calculated
by
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md),
if the feeding level is nonzero this can also be expressed as
\$\$E\_{r.i}(w) = \alpha_i\\ {\tt feeding\\level}\_i(w)\\ h_i(w) - {\tt
metab}\_i(w)\$\$ where \\h_i\\ is the maximum intake rate, set with
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md).
However this function is using the first equation above so that it works
also when the maximum intake rate is infinite, i.e., there is no
satiation. The assimilation rate \\\alpha_i\\ is taken from the species
parameter data frame in `params`. The metabolic rate `metab` is taken
from `params` and set with
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md).

The return value can be negative, which means that the energy intake
does not cover the cost of metabolism and movement.

## Your own energy rate function

By default
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md)
calls `mizerEReproAndGrowth()`. However you can replace this with your
own alternative energy rate function. If your function is called
`"myEReproAndGrowth"` then you register it in a MizerParams object
`params` with

    params <- setRateFunction(params, "EReproAndGrowth", "myEReproAndGrowth")

Your function will then be called instead of `mizerEReproAndGrowth()`,
with the same arguments.

## See also

Other mizer rate functions:
[`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md),
[`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md),
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
