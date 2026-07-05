# Set up parameters for a single species in a power-law background

**\[experimental\]**

This functions creates a `MizerParams` object with a single species.
This species is embedded in a fixed power-law community spectrum
\$\$N_c(w) = \kappa w^{-\lambda}\$\$ This community provides the food
income for the species. Cannibalism is switched off. The predation
mortality arises only from the predators in the power-law community and
it is assumed that the predators in the community have the same feeding
parameters as the foreground species. The function has many arguments,
all of which have default values.

## Usage

``` r
newSingleSpeciesParams(
  species_name = "Target species",
  w_max = 100,
  w_min = 0.001,
  eta = 10^(-0.6),
  w_mat = w_max * eta,
  no_w = log10(w_max/w_min) * 20 + 1,
  n = 3/4,
  p = n,
  lambda = 2.05,
  kappa = 0.005,
  alpha = 0.4,
  h = 30,
  beta = 100,
  sigma = 1.3,
  f0 = 0.6,
  fc = 0.25,
  ks = NA,
  gamma = NA,
  ext_mort_prop = 0,
  reproduction_level = 0,
  second_order_w = FALSE,
  R_factor = deprecated(),
  w_inf = deprecated(),
  k_vb = deprecated()
)
```

## Arguments

- species_name:

  A string with a name for the species. Will be used in plot legends.

- w_max:

  Maximum size of species

- w_min:

  Egg size of species

- eta:

  Ratio between maturity size `w_mat` and maximum size `w_max`. Default
  is 10^(-0.6), approximately 1/4. Ignored if `w_mat` is supplied
  explicitly.

- w_mat:

  Maturity size of species. Default value is `eta * w_max`.

- no_w:

  The number of size bins in the community spectrum. These bins will be
  equally spaced on a logarithmic scale. Default value is such that
  there are 20 bins for each factor of 10 in weight.

- n:

  Scaling exponent of the maximum intake rate.

- p:

  Scaling exponent of the standard metabolic rate. By default this is
  equal to the exponent `n`.

- lambda:

  Exponent of the abundance power law.

- kappa:

  Coefficient in abundance power law.

- alpha:

  The assimilation efficiency.

- h:

  Maximum food intake rate.

- beta:

  Preferred predator prey mass ratio.

- sigma:

  Width of prey size preference.

- f0:

  Expected average feeding level. Used to set `gamma`, the coefficient
  in the search rate. Ignored if `gamma` is given explicitly.

- fc:

  Critical feeding level. Used to determine `ks` if it is not given
  explicitly.

- ks:

  Standard metabolism coefficient. If not provided, default will be
  calculated from critical feeding level argument `fc`.

- gamma:

  Volumetric search rate. If not provided, default is determined by
  [`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md)
  using the value of `f0`.

- ext_mort_prop:

  The proportion of the total mortality that comes from external
  mortality, i.e., from sources not explicitly modelled. A number in the
  interval \[0, 1).

- reproduction_level:

  A number between 0 and 1 that determines the level of density
  dependence in reproduction, see
  [`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md).

- second_order_w:

  **\[experimental\]** Selects the second-order numerical scheme for the
  new model, applied before the resource and abundance power laws are
  constructed so that they are built bin-averaged from the start.
  Accepts the same values as the
  [`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
  setter. Defaults to `FALSE`. See
  [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md).

- R_factor:

  **\[deprecated\]** Use `reproduction_level = 1 / R_factor` instead.

- w_inf:

  **\[deprecated\]** The argument has been renamed to `w_max`.

- k_vb:

  **\[deprecated\]** The von Bertalanffy growth parameter.

## Value

An object of type `MizerParams`

## Details

In addition to setting up the parameters, this function also sets up an
initial condition that is close to steady state, under the assumption of
no fishing.

The function rounds `no_w` to the nearest integer and increases it if
necessary so that there are at least 5 size bins per factor 10 in body
size. It requires `w_min < w_mat < w_max`, `ext_mort_prop` in `[0, 1)`,
positive values for `n`, `lambda`, `kappa`, `alpha`, `h`, `beta`,
`sigma` and `f0`, and `fc` between 0 and `f0` if `fc` is supplied. If
`gamma` is supplied then `f0` is ignored. The function stops if the
resulting feeding level is not sufficient to maintain the species.

The returned model has a single foreground species with cannibalism
switched off and a fixed power-law background community that provides
both food and predation mortality. The initial species spectrum is
scaled so that its maximum abundance is half the background abundance at
the corresponding size, and `erepro` is then adjusted so the initial
state satisfies the egg boundary condition.

The diffusion rate is set to `0`

## See also

Other functions for setting up models:
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md),
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md),
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)

## Examples

``` r
params <- newSingleSpeciesParams()
sim <- project(params, t_max = 5, effort = 0)
plotSpectra(sim)
```
