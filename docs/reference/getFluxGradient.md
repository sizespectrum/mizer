# Get flux gradient

Calculates the flux divergence \\(J\_{j+1} - J_j)/\Delta w_j\\ that
appears as the second term in the discretised size-spectrum transport
equation \$\$\frac{\partial N_j}{\partial t} + \frac{J\_{j+1} -
J_j}{\Delta w_j} = -\mu_j N_j.\$\$ The bin-boundary fluxes \\J_j\\ are
obtained from
[`getFlux()`](https://sizespectrum.org/mizer/reference/getFlux.md),
which uses the advective-flux scheme stored in the `flux` entry of the
[`second_order_w`](https://sizespectrum.org/mizer/reference/second_order_w.md)
slot of `params`. The flux leaving the largest size class through the
upper boundary (\\J\_{K+1}\\) is evaluated with the same scheme using
the boundary condition \\N\_{K+1} = 0\\.

## Usage

``` r
getFluxGradient(object, ...)
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

- `MizerParams`: An `ArraySpeciesBySize` object (species x size) giving
  the flux divergence in each size bin, in units of \\g^{-1} \\
  \text{year}^{-1}\\.

- `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species
  x size) with the flux divergence at every saved time step. If
  `drop = TRUE` then dimensions of length 1 will be removed.

## See also

[`getFlux()`](https://sizespectrum.org/mizer/reference/getFlux.md),
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)

Other rate functions:
[`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md),
[`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md),
[`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md),
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md),
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md),
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md),
[`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md),
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md),
[`getFlux()`](https://sizespectrum.org/mizer/reference/getFlux.md),
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md),
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md),
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md),
[`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md),
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md),
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md),
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

## Examples

``` r
# \donttest{
params <- NS_params
fg <- getFluxGradient(params)
sim <- project(params, t_max = 5)
fg_sim <- getFluxGradient(sim)
# }
```
