# Example species parameter set based on the North Sea

This data set is based on species in the North Sea (Blanchard et al.).
It is a data.frame that contains all the necessary information to be
used by the
[`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.md)
constructor. As there is no gear column, each species is assumed to be
fished by a separate gear.

## Usage

``` r
NS_species_params
```

## Format

A data frame with 12 rows and 7 columns. Each row is a species.

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

## Source

Blanchard et al.

## Examples

``` r
params <- newMultispeciesParams(NS_species_params)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No h provided for some species, so using age at maturity to calculate it.
#> ℹ Because you have n != p, the default value for `h` is not very good.
#> ℹ Because the age at maturity is not known, I need to fall back to using
#> von Bertalanffy parameters, where available.
#> ℹ Using z0 = z0pre * w_inf ^ z0exp for calculated z0 values.
#> ℹ Using f0, h, lambda, kappa and the predation kernel to calculate gamma.
```
