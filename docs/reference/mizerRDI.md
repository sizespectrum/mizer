# Get density-independent rate of reproduction needed to project standard mizer model

Calculates the density-independent rate of total egg production
\\R\_{di}\\ (units 1/year) before density dependence, by species. You
would not usually call this function directly but instead use
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md), which
then calls this function unless an alternative function has been
registered, see below.

## Usage

``` r
projectRDI(
  params,
  n,
  n_pp,
  n_other,
  t = 0,
  e_growth,
  mort,
  e_repro,
  diffusion = NULL,
  ...
)

# S3 method for class 'MizerParams'
projectRDI(
  params,
  n,
  n_pp,
  n_other,
  t = 0,
  e_growth,
  mort,
  e_repro,
  diffusion = NULL,
  ...
)

mizerRDI(
  params,
  n,
  n_pp,
  n_other,
  t = 0,
  e_growth,
  mort,
  e_repro,
  diffusion = NULL,
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

- e_growth:

  An array (species x size) with the energy available for growth as
  calculated by
  [`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md).
  Unused.

- mort:

  An array (species x size) with the mortality rate as calculated by
  [`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md).
  Unused.

- e_repro:

  An array (species x size) with the energy available for reproduction
  as calculated by
  [`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md).

- diffusion:

  An array (species x size) with the diffusion rate as calculated by
  [`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md).
  Unused by the default function but supplied to custom RDI functions.

- ...:

  Unused

## Value

A numeric vector with the rate of egg production for each species.

## Details

This rate is obtained by taking the per capita rate \\E_r(w)\psi(w)\\ at
which energy is invested in reproduction, as calculated by
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
multiplying it by the number of individuals\\N(w)\\ and integrating over
all sizes \\w\\ and then multiplying by the reproductive efficiency
\\\epsilon\\ and dividing by the egg size `w_min`, and by a factor of
two to account for the two sexes: \$\$R\_{di} = \frac{\epsilon}{2
w\_{min}} \int N(w) E_r(w) \psi(w) \\ dw\$\$

Used by [`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md)
to calculate the actual, density dependent rate. See
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
for more details.

## Your own reproduction function

By default
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md) calls
`mizerRDI()`. However you can replace this with your own alternative
reproduction function. If your function is called `"myRDI"` then you
register it in a MizerParams object `params` with

    params <- setRateFunction(params, "RDI", "myRDI")

Your function will then be called instead of `mizerRDI()`, with the same
arguments. For an example of an alternative reproduction function see
[`constantEggRDI()`](https://sizespectrum.org/mizer/reference/constantEggRDI.md).

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
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md),
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
