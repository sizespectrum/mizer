# Get encounter rate

Returns the rate at which a predator of species \\i\\ and weight \\w\\
encounters food (grams/year).

## Usage

``` r
getEncounter(object, ...)
```

## Arguments

- object:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  or
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  object.

- ...:

  Additional arguments that depend on the class of `object`.

  **For a
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object:**

  `n`

  :   A matrix of species abundances (species x size). Defaults to the
      initial abundances stored in `object`.

  `n_pp`

  :   A vector of the resource abundance by size. Defaults to the
      initial resource abundance stored in `object`.

  `n_other`

  :   A named list of the abundances of other dynamical components.
      Defaults to the initial values stored in `object`.

  `t`

  :   The time for which to do the calculation. Defaults to 0.

  **For a
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  object:**

  `time_range`

  :   The time range over which to return the rates. Either a vector of
      values, a vector of min and max time, or a single value. Defaults
      to the whole time range of the simulation.

  `drop`

  :   If `TRUE` then any dimension of length 1 is removed from the
      returned array.

## Value

- `MizerParams`: An `ArraySpeciesBySize` object (predator species x
  predator size) with the encounter rates.

- `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x predator
  species x predator size) with the encounter rates at every time step.
  If `drop = TRUE` then dimensions of length 1 will be removed.

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

[`projectEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
is the S3 generic used by extension-aware projections. Extension
packages can add methods for their marker classes and call
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) to compose
encounter-rate changes. The `MizerParams` method contains the standard
mizer calculation and is also exported as
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
for compatibility.

## Your own encounter function

By default `getEncounter()` calls
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
on models without extensions. However you can replace this with your own
alternative encounter function. If your function is called
`"myEncounter"` then you register it in a MizerParams object `params`
with

    params <- setRateFunction(params, "Encounter", "myEncounter")

Your function will then be called instead of
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md),
with the same arguments.

## See also

Other rate functions:
[`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md),
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md),
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md),
[`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md),
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md),
[`getFlux()`](https://sizespectrum.org/mizer/reference/getFlux.md),
[`getFluxGradient()`](https://sizespectrum.org/mizer/reference/getFluxGradient.md),
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md),
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md),
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md),
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md),
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md),
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
encounter <- getEncounter(NS_params)
str(encounter)
#>  'ArraySpeciesBySize' num [1:12, 1:100] 0.299 0.453 0.502 0.575 0.492 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ sp: chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>  - attr(*, "value_name")= chr "Encounter rate"
#>  - attr(*, "units")= chr "g/year"
#>  - attr(*, "representation")= chr "point"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
```
