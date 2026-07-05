# Resource parameters

The recommended way to change the resource dynamics parameters is to use
[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md).
The `resource_params` list contains values that are helpful in setting
up the actual size-dependent parameters with
[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md).
If you have specified a custom resource dynamics function that requires
additional parameters, then these should also be added to the
`resource_params` list.

## Usage

``` r
resource_params(params)

resource_params(params) <- value
```

## Arguments

- params:

  A MizerParams object

- value:

  A named list of resource parameters.

## Value

A named list of resource parameters.

## Details

The `resource_params` list will at least contain the slots `kappa`,
`lambda`, `w_pp_cutoff` and `n`.

The resource parameter `n` is the exponent for the power-law form for
the replenishment rate \\r_R(w)\\: \$\$r_R(w) = r_R\\ w^{n-1}.\$\$

The resource parameter `lambda` (\\\lambda\\) is the exponent for the
power-law form for the carrying capacity \\c_R(w)\\ and `w_pp_cutoff` is
its cutoff value: \$\$c_R(w) = c_R w^{-\lambda}\$\$ for all \\w\\ less
than `w_pp_cutoff` and zero for larger sizes.

The resource parameter `kappa` (\\\kappa\\) is slightly different in
that it is not a parameter for the resource dynamics. Instead it is a
parameter that determined the initial resource abundance when the model
was created: \$\$N_R(w) = \kappa\\ w^{-\lambda}\$\$ for all \\w\\ less
than `w_pp_cutoff` and zero for larger sizes. Note that the initial
resource abundance is not changed by
[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md)
even if you change the value of `kappa` in the `resource_params`.

## See also

[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md)
