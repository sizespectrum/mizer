# Deprecated function for setting up parameters for a trait-based model

**\[deprecated\]**

This function has been deprecated in favour of the function
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)
that sets better default values.

## Usage

``` r
set_trait_model(
  no_sp = 10,
  min_w_inf = 10,
  max_w_inf = 1e+05,
  no_w = 100,
  min_w = 0.001,
  max_w = max_w_inf * 1.1,
  min_w_pp = 1e-10,
  w_pp_cutoff = 1,
  k0 = 50,
  n = 2/3,
  p = 0.75,
  q = 0.9,
  eta = 0.25,
  r_pp = 4,
  kappa = 0.005,
  lambda = 2 + q - n,
  alpha = 0.6,
  ks = 4,
  z0pre = 0.6,
  h = 30,
  beta = 100,
  sigma = 1.3,
  f0 = 0.5,
  gamma = NA,
  knife_edge_size = 1000,
  gear_names = "knife_edge_gear",
  ...
)
```

## Arguments

- no_sp:

  The number of species in the model. The default value is 10. The more
  species, the longer takes to run.

- min_w_inf:

  The asymptotic size of the smallest species in the community.

- max_w_inf:

  The asymptotic size of the largest species in the community.

- no_w:

  The number of size bins in the community spectrum.

- min_w:

  The smallest size of the community spectrum.

- max_w:

  Maximum size of the consumer size grid passed to
  [`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.md).
  Default value is `max_w_inf * 1.1`.

- min_w_pp:

  Smallest size on the resource size grid passed to
  [`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.md).
  Default value is `1e-10`.

- w_pp_cutoff:

  The cut off size of the resource spectrum. Default value is 1.

- k0:

  Multiplier for the maximum recruitment. Default value is 50.

- n:

  Scaling of the intake. Default value is 2/3.

- p:

  Scaling of the standard metabolism. Default value is 0.75.

- q:

  Exponent of the search volume. Default value is 0.9.

- eta:

  Factor to calculate `w_mat` from asymptotic size.

- r_pp:

  Growth rate parameter for the resource spectrum. Default value is 4.

- kappa:

  Coefficient in abundance power law. Default value is 0.005.

- lambda:

  Exponent of the abundance power law. Default value is (2+q-n).

- alpha:

  The assimilation efficiency of the community. The default value is 0.6

- ks:

  Standard metabolism coefficient. Default value is 4.

- z0pre:

  The coefficient of the background mortality of the community. z0 =
  z0pre \* w_inf ^ (n-1). The default value is 0.6.

- h:

  Maximum food intake rate. Default value is 30.

- beta:

  Preferred predator prey mass ratio. Default value is 100.

- sigma:

  Width of prey size preference. Default value is 1.3.

- f0:

  Expected average feeding level. Used to set `gamma`, the factor for
  the search volume. The default value is 0.5.

- gamma:

  Volumetric search rate. Estimated using `h`, `f0` and `kappa` if not
  supplied.

- knife_edge_size:

  The minimum size at which the gear or gears select species. Must be of
  length 1 or no_sp.

- gear_names:

  The names of the fishing gears. A character vector, the same length as
  the number of species. Default is 1 - no_sp.

- ...:

  Other arguments to pass to the `MizerParams` constructor.

## Value

An object of type `MizerParams`

## Details

This functions creates a `MizerParams` object so that trait-based-type
models can be easily set up and run. The trait-based size spectrum model
can be derived as a simplification of the general size-based model used
in `mizer`. The species-specific parameters are the same for all
species, except for the asymptotic size, which is considered the most
important trait characterizing a species. Other parameters are related
to the asymptotic size. For example, the size at maturity is given by
`w_max * eta`, where `eta` is the same for all species. For the
trait-based model the number of species is not important. For
applications of the trait-based model see Andersen & Pedersen (2010).
See the `mizer` vignette for more details and examples of the
trait-based model.

The function has many arguments, all of which have default values. Of
particular interest to the user are the number of species in the model
and the minimum and maximum asymptotic sizes. The asymptotic sizes of
the species are spread evenly on a logarithmic scale within this range.

The stock recruitment relationship is the default Beverton-Holt style.
The maximum recruitment is calculated using equilibrium theory (see
Andersen & Pedersen, 2010) and a multiplier, `k0`. Users should adjust
`k0` to get the spectra they want.

The factor for the search volume, `gamma`, is calculated using the
expected feeding level, `f0`.

Fishing selectivity is modelled as a knife-edge function with one
parameter, `knife_edge_size`, which is the size at which species are
selected. Each species can either be fished by the same gear
(`knife_edge_size` has a length of 1) or by a different gear (the length
of `knife_edge_size` has the same length as the number of species and
the order of selectivity size is that of the asymptotic size).

The resulting `MizerParams` object can be projected forward using
`project` like any other `MizerParams` object. When projecting the
community model it may be necessary to reduce `dt` to 0.1 to avoid any
instabilities with the solver. You can check this by plotting the
biomass or abundance through time after the projection.

## References

K. H. Andersen and M. Pedersen, 2010, Damped trophic cascades driven by
fishing in model marine ecosystems. Proceedings of the Royal Society V,
Biological Sciences, 1682, 795-802.
