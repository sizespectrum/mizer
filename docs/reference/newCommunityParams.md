# Set up parameters for a community-type model

This functions creates a
[`MizerParams`](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
object describing a community-type model. The function has many
arguments, all of which have default values.

## Usage

``` r
newCommunityParams(
  max_w = 1e+06,
  min_w = 0.001,
  no_w = 100,
  min_w_pp = 1e-10,
  z0 = 0.1,
  alpha = 0.2,
  f0 = 0.7,
  h = 10,
  gamma = NA,
  beta = 100,
  sigma = 2,
  n = 2/3,
  kappa = 1000,
  lambda = 2.05,
  r_pp = 10,
  knife_edge_size = 1000,
  reproduction,
  second_order_w = FALSE,
  info_level = 2
)
```

## Arguments

- max_w:

  The maximum size of the community. The `w_max` of the species used to
  represent the community is set to this value.

- min_w:

  The minimum size of the community.

- no_w:

  The number of size bins in the consumer spectrum.

- min_w_pp:

  The smallest size of the resource spectrum. By default this is set to
  the smallest value at which any of the consumers can feed.

- z0:

  The background mortality of the community.

- alpha:

  The assimilation efficiency of the community.

- f0:

  The average feeding level of individuals who feed on a power-law
  spectrum. This value is used to calculate the search rate parameter
  `gamma`.

- h:

  The coefficient of the maximum food intake rate.

- gamma:

  Volumetric search rate. Passed through to
  [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md),
  which estimates it from `h`, `f0`, and `kappa` if it is left as `NA`.

- beta:

  The preferred predator prey mass ratio.

- sigma:

  The width of the prey preference.

- n:

  The allometric growth exponent. Used as allometric exponent for the
  maximum intake rate of the community as well as the intrinsic growth
  rate of the resource.

- kappa:

  The coefficient of the initial resource abundance power-law.

- lambda:

  Used to set power-law exponent for resource capacity if the
  `resource_capacity` argument is given as a single number.

- r_pp:

  Growth rate parameter for the resource spectrum.

- knife_edge_size:

  The size at the edge of the knife-edge-selectivity function.

- reproduction:

  The constant reproduction in the smallest size class of the community
  spectrum. By default this is set to the rate required to maintain the
  constructed initial egg abundance.

- second_order_w:

  **\[experimental\]** Selects the second-order numerical scheme for the
  new model. Accepts the same values as the
  [`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
  setter: a single logical (`TRUE` switches on both second-order flux
  and bin-averaging), a single flux scheme name (`"upwind"`,
  `"van_leer"` or `"centred"`), or a named vector with entries `flux`
  and/or `bin_average`. The `bin_average` choice is applied *before* the
  resource and abundance power laws are constructed, so they are built
  bin-averaged from the start (unlike setting
  [`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
  on an existing object). The `flux` scheme governs time projection
  only, so the robust first-order upwind scheme is used for the
  construction-time steady-state solve and the chosen scheme is then
  activated for the returned model. Defaults to `FALSE` (the first-order
  behaviour of previous mizer).

- info_level:

  Controls the amount of information messages that are shown when the
  function sets default values for parameters. Higher levels lead to
  more messages.

## Value

An object of type
[`MizerParams`](https://sizespectrum.org/mizer/reference/MizerParams-class.md)

## Details

A community model has several features that distinguish it from a
multi-species model:

- Species identities of individuals are ignored. All are aggregated into
  a single community.

- The resource spectrum only extends to the start of the community
  spectrum.

- Reproductive rate is constant, independent of the energy invested in
  reproduction, which is set to 0.

- Standard metabolism is turned off (the parameter `ks` is set to 0).
  Consequently, the growth rate is now determined solely by the
  assimilated food

Fishing selectivity is modelled as a knife-edge function with one
parameter, `knife_edge_size`, which determines the size at which species
are selected.

The resulting `MizerParams` object can be projected forward using
[`project()`](https://sizespectrum.org/mizer/reference/project.md) like
any other `MizerParams` object. When projecting the community model it
may be necessary to keep a small time step size `dt` of around 0.1 to
avoid any instabilities with the solver. You can check for these
numerical instabilities by plotting the biomass or abundance through
time after the projection.

## References

K. H. Andersen,J. E. Beyer and P. Lundberg, 2009, Trophic and individual
efficiencies of size-structured communities, Proceedings of the Royal
Society, 276, 109-114

## See also

Other functions for setting up models:
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md),
[`newSingleSpeciesParams()`](https://sizespectrum.org/mizer/reference/newSingleSpeciesParams.md),
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)

## Examples

``` r
params <- newCommunityParams()
sim <- project(params, t_max = 10)
plotBiomass(sim)

plotSpectra(sim, power = 2)


# More satiation. More mortality
params <- newCommunityParams(f0 = 0.8, z0 = 0.4)
sim <- project(params, t_max = 10)
plotBiomass(sim)

plotSpectra(sim, power = 2)
```
