# Get diffusion rate from predation

Calculates the diffusion rate \\D_i(w)\\ (grams^2/year) for each
species. This diffusion rate has two components:

1.  The diffusion due due to the variability in prey sizes. This is the
    diffusion term from the jump-growth equation.

2.  Any externally specified diffusion, which is added via
    [`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md)

## Usage

``` r
getDiffusion(object, ...)
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
  predator size) with the diffusion rates.

- `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x predator
  species x predator size) with the diffusion rates at every time step.
  If `drop = TRUE` then dimensions of length 1 will be removed.

## Details

The diffusion due due to the variability in prey sizes is determined by
summing over all prey species and the resource spectrum and then
integrating over all prey sizes \\w_p\\, weighted by predation kernel
\\\phi(w,w_p)\\: \$\$ d_i(w) =
(1-f_i(w))(\alpha_i(1-\psi_i(w)))^2\gamma_i(w) \int \left( \theta\_{ip}
N_R(w_p) + \sum\_{j} \theta\_{ij} N_j(w_p) \right) \phi_i(w,w_p) w_p^2
\\ dw_p. \$\$ Here \\N_j(w)\\ is the abundance density of species \\j\\
and \\N_R(w)\\ is the abundance density of resource. The overall
prefactor \\\gamma_i(w)\\ determines the predation power of the
predator. It could be interpreted as a search volume and is set with the
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md)
function. The predation kernel \\\phi(w,w_p)\\ is set with the
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
function. The species interaction matrix \\\theta\_{ij}\\ is set with
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md)
and the resource interaction vector \\\theta\_{ip}\\ is taken from the
`interaction_resource` column in
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md).
\\f(w)\\ is the feeding level calculated with
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md).
\\\psi(w)\\ is the proportion of the available energy that is invested
in reproduction instead of growth, obtained with
[`psi()`](https://sizespectrum.org/mizer/reference/setReproduction.md).

The diffusion integral is normally evaluated efficiently with a fast
Fourier transform, which assumes that the predation kernel depends only
on the ratio of predator to prey size. If a custom predation kernel that
depends on predator and prey size separately has been set with
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
the integral is instead evaluated by direct summation over the full
predation kernel, as in
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md).

## References

Datta, S., Delius, G. W. and Law, R. (2010). A jump-growth model for
predator-prey dynamics: derivation and application to marine ecosystems.
Bulletin of Mathematical Biology, 72(6):1361–1382

## See also

Other rate functions:
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md),
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md),
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
