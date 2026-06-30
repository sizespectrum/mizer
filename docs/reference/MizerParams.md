# Alias for `set_multispecies_model()`

**\[deprecated\]** An alias provided for backward compatibility with
mizer version \<= 1.0

## Usage

``` r
MizerParams(
  species_params,
  interaction = matrix(1, nrow = nrow(species_params), ncol = nrow(species_params)),
  min_w_pp = 1e-10,
  min_w = 0.001,
  max_w = NULL,
  no_w = 100,
  n = 2/3,
  q = 0.8,
  f0 = 0.6,
  kappa = 1e+11,
  lambda = 2 + q - n,
  r_pp = 10,
  ...
)
```

## Arguments

- species_params:

  A data frame of species-specific parameter values.

- interaction:

  Optional interaction matrix of the species (predator species x prey
  species). By default all entries are 1. See "Setting interaction
  matrix" section below.

- min_w_pp:

  The smallest size of the resource spectrum. By default this is set to
  the smallest value at which any of the consumers can feed.

- min_w:

  Sets the size of the eggs of all species for which this is not given
  in the `w_min` column of the `species_params` dataframe.

- max_w:

  The largest size of the consumer spectrum. By default this is set to
  the largest `w_max` specified in the `species_params` data frame.

- no_w:

  The number of size bins in the consumer spectrum.

- n:

  The allometric growth exponent. This can be overruled for individual
  species by including a `n` column in the `species_params`.

- q:

  Allometric exponent of search volume

- f0:

  Expected average feeding level. Used to set `gamma`, the coefficient
  in the search rate. Ignored if `gamma` is given explicitly.

- kappa:

  The coefficient of the initial resource abundance power-law.

- lambda:

  Used to set power-law exponent for resource capacity if the
  `resource_capacity` argument is given as a single number.

- r_pp:

  **\[deprecated\]**. Use `resource_rate` argument instead.

- ...:

  Further arguments passed to
  [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md).

## Value

A MizerParams object

## Details

If `species_params` contains a `w_inf` column then it is copied to
`w_max`. If `max_w` is not supplied then it is set to
`1.1 * max(species_params$w_max)`. The supplied `min_w_pp` is shifted up
by one grid step before being passed to
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
to compensate for the fact that newer mizer versions extend the full
size grid below `min_w_pp`.

Missing legacy columns in `species_params` are filled as follows:
`gear = species`, `k = 0`, `alpha = 0.6`, `erepro = 1`,
`sel_func = "knife_edge"`, `knife_edge_size = w_mat` if needed,
`catchability = 1`, `ks = h * 0.2`, and `m = 1`. If `h` is missing it is
calculated from `k_vb`, `alpha`, `f0` and `w_max`. If `gamma` is missing
it is calculated from `f0`, `h`, `beta`, `sigma`, `lambda` and `kappa`.
