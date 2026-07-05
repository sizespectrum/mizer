# Calculate age at maturity

Uses the size-dependent growth rate and the size at maturity to
calculate the age at maturity.

## Usage

``` r
age_mat(params, ...)
```

## Arguments

- params:

  A MizerParams object

- ...:

  Currently unused.

## Value

A named vector. The names are the species names and the values are the
ages at maturity.

## Details

Using that by definition of the growth rate \\g(w) = dw/dt\\ we have
that \$\$\mathrm{age\_{mat}} = \int_0^{w\_{mat}.}\frac{dw}{g(w)}\$\$ In
the implementation this integral is approximated on the model size grid
by summing `dw / g(w)` over all size bins with `w < w_mat`.

## Examples

``` r
age_mat(NS_params)
#>     Sprat   Sandeel    N.pout   Herring       Dab   Whiting      Sole   Gurnard 
#> 2.5560864 0.9832776 1.7962884 3.0382630 1.9612689 2.7122876 4.4940829 4.1384558 
#>    Plaice   Haddock       Cod    Saithe 
#> 6.2608150 3.4391147 2.8583434 4.6315535 
```
