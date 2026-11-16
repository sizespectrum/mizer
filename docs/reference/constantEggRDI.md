# Choose egg production to keep egg density constant

**\[experimental\]** The new egg production is set to compensate for the
loss of individuals from the smallest size class through growth and
mortality. The result should not be modified by density dependence, so
this should be used together with the
[`noRDD()`](https://sizespectrum.org/mizer/reference/noRDD.md) function,
see example.

## Usage

``` r
constantEggRDI(params, n, e_growth, mort, diffusion, ...)
```

## Arguments

- params:

  A MizerParams object

- n:

  A matrix of species abundances (species x size).

- e_growth:

  A two dimensional array (species x size) holding the energy available
  for growth as calculated by
  [`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md).

- mort:

  A two dimensional array (species x size) holding the mortality rate as
  calculated by
  [`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md).

- diffusion:

  A two dimensional array (species x size) holding the diffusion rate as
  calculated by
  [`mizerDiffusion()`](https://sizespectrum.org/mizer/reference/mizerDiffusion.md).

- ...:

  Unused

## Value

Vector with the value for each species

## See also

Other functions calculating density-dependent reproduction rate:
[`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md),
[`RickerRDD()`](https://sizespectrum.org/mizer/reference/RickerRDD.md),
[`SheperdRDD()`](https://sizespectrum.org/mizer/reference/SheperdRDD.md),
[`constantRDD()`](https://sizespectrum.org/mizer/reference/constantRDD.md),
[`noRDD()`](https://sizespectrum.org/mizer/reference/noRDD.md)

## Examples

``` r
# \donttest{
# choose an example params object
params <- NS_params
# We set the reproduction rate functions
params <- setRateFunction(params, "RDI", "constantEggRDI")
params <- setRateFunction(params, "RDD", "noRDD")
# Now the egg density should stay fixed no matter how we fish
sim <- project(params, effort = 10, progress_bar = FALSE)
# To check that indeed the egg densities have not changed, we first construct
# the indices for addressing the egg densities
no_sp <- nrow(params@species_params)
idx <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
# Now we can check equality between egg densities at the start and the end
all.equal(finalN(sim)[idx], initialN(params)[idx])
#> [1] TRUE
# }
```
