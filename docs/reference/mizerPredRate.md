# Get predation rate needed to project standard mizer model

Calculates the potential rate (in units 1/year) at which a prey
individual of a given size \\w\\ is killed by predators from species
\\j\\. In formulas \$\${\tt pred\\rate}\_j(w_p) = \int \phi_j(w,w_p)
(1-f_j(w)) \gamma_j(w) N_j(w) \\ dw.\$\$ This potential rate is used in
the function
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md)
to calculate the realised predation mortality rate on the prey
individual. You would not usually call this function directly but
instead use
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md),
which then calls this function unless an alternative function has been
registered, see below.

## Usage

``` r
projectPredRate(params, n, n_pp, n_other, t = 0, feeding_level, ...)

# S3 method for class 'MizerParams'
projectPredRate(params, n, n_pp, n_other, t = 0, feeding_level, ...)

mizerPredRate(params, n, n_pp, n_other, t = 0, feeding_level, ...)
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

- feeding_level:

  An array (species x size) with the feeding level as calculated by
  [`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md).

- ...:

  Unused

## Value

A named two dimensional array (predator species x prey size) with the
predation rate, where the prey size runs over fish community plus
resource spectrum.

## Your own predation rate function

By default
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md)
calls `mizerPredRate()`. However you can replace this with your own
alternative predation rate function. If your function is called
`"myPredRate"` then you register it in a MizerParams object `params`
with

    params <- setRateFunction(params, "PredRate", "myPredRate")

Your function will then be called instead of `mizerPredRate()`, with the
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
[`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md),
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md),
[`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md),
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md),
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
