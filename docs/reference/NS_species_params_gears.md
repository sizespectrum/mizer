# Example species parameter set based on the North Sea with different gears

This data set is based on species in the North Sea (Blanchard et al.).
It is similar to the data set `NS_species_params` except that this one
has an additional column specifying the fishing gear that operates on
each species.

## Usage

``` r
NS_species_params_gears
```

## Format

A data frame with 12 rows and 8 columns. Each row is a species.

- species:

  Name of the species

- w_max:

  Computational upper size boundary, defaulting to `1.5 * w_inf`.

- w_mat:

  Size at maturity

- beta:

  Size preference ratio

- sigma:

  Width of the size-preference

- R_max:

  Maximum reproduction rate

- k_vb:

  The von Bertalanffy k parameter

- w_inf:

  The von Bertalanffy asymptotic size. This is the required maximum-size
  parameter.

- gear:

  Name of the fishing gear

## Source

Blanchard et al.

## Examples

``` r
params <- MizerParams(NS_species_params_gears)
#> Warning: `set_multispecies_model()` was deprecated in mizer 2.0.0.
#> ℹ Please use `newMultispeciesParams()` instead.
#> Note: No sel_func column in species data frame. Setting selectivity to be 'knife_edge' for all species.
#> Note: No knife_edge_size column in species data frame. Setting knife edge selectivity equal to w_mat.
#> Note: No h column in species data frame so using f0 and k_vb to calculate it.
#> Note: No gamma column in species data frame so using f0, h, beta, sigma, lambda and kappa to calculate it.
#> Note: No ks column in species data frame so using ks = h * 0.2.
#> Note: No m column in species data frame so using m = 1.
#> Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.
```
