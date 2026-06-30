# Example interaction matrix for the North Sea example

The interaction coefficient between predator and prey species in the
North Sea.

## Usage

``` r
NS_interaction
```

## Format

A 12 x 12 matrix.

## Source

Blanchard et al.

## Examples

``` r
params <- newMultispeciesParams(NS_species_params_gears,
                                interaction = NS_interaction)
#> Because you have n != p, the default value for `h` is not very good.
#> Because the age at maturity is not known, I need to fall back to using
#> von Bertalanffy parameters, where available, and this is not reliable.
#> No ks column so calculating from critical feeding level.
#> Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.
#> Using f0, h, lambda, kappa and the predation kernel to calculate gamma.
```
