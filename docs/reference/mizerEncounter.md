# Get encounter rate during projection

Calculates the rate \\E_i(w)\\ at which a predator of species \\i\\ and
weight \\w\\ encounters food (grams/year). You would not usually call
this function directly but instead use
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md),
which then calls this function unless an alternative function has been
registered, see below.

## Usage

``` r
projectEncounter(params, n, n_pp, n_other, t = 0, ...)

# S3 method for class 'MizerParams'
projectEncounter(params, n, n_pp, n_other, t = 0, ...)

mizerEncounter(params, n, n_pp, n_other, t = 0, ...)
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

- ...:

  Unused

## Value

A named two dimensional array (predator species x predator size) with
the encounter rates.

## Predation encounter

The encounter rate \\E_i(w)\\ at which a predator of species \\i\\ and
weight \\w\\ encounters food has contributions from the encounter of
fish prey and of resource. This is determined by summing over all prey
species and the resource spectrum and then integrating over all prey
sizes \\w_p\\, weighted by predation kernel \\\phi(w,w_p)\\: \$\$ E_i(w)
= \gamma_i(w) \int \left( \theta\_{ip} N_R(w_p) + \sum\_{j} \theta\_{ij}
N_j(w_p) \right) \phi_i(w,w_p) w_p \\ dw_p. \$\$ Here \\N_j(w)\\ is the
abundance density of species \\j\\ and \\N_R(w)\\ is the abundance
density of resource. The overall prefactor \\\gamma_i(w)\\ determines
the predation power of the predator. It could be interpreted as a search
volume and is set with the
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md)
function. The predation kernel \\\phi(w,w_p)\\ is set with the
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
function. The species interaction matrix \\\theta\_{ij}\\ is set with
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md)
and the resource interaction vector \\\theta\_{ip}\\ is taken from the
`interaction_resource` column in `params@species_params`.

## Details

The encounter rate is multiplied by \\1-f_0\\ to obtain the consumption
rate, where \\f_0\\ is the feeding level calculated with
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md).
This is used by the
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
function for performing simulations.

The function returns values also for sizes outside the size-range of the
species. These values should not be used, as they are meaningless.

If your model contains additional components that you added with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
and for which you specified an `encounter_fun` function then the
encounters of these components will be included in the returned value.

## Extension hook

`projectEncounter()` is the S3 generic used by extension-aware
projections. Extension packages can add methods for their marker classes
and call [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) to
compose encounter-rate changes. The `MizerParams` method contains the
standard mizer calculation and is also exported as `mizerEncounter()`
for compatibility.

## Your own encounter function

By default
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
calls `mizerEncounter()` on models without extensions. However you can
replace this with your own alternative encounter function. If your
function is called `"myEncounter"` then you register it in a MizerParams
object `params` with

    params <- setRateFunction(params, "Encounter", "myEncounter")

Your function will then be called instead of `mizerEncounter()`, with
the same arguments.

## See also

Other mizer rate functions:
[`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md),
[`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md),
[`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md),
[`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md),
[`mizerFMortGear()`](https://sizespectrum.org/mizer/reference/mizerFMortGear.md),
[`mizerFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md),
[`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md),
[`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md),
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md),
[`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md),
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md),
[`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
