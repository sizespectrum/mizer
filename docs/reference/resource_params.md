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

The resource parameter `kappa` (\\\kappa\\) is the coefficient \\c_R\\
of the carrying capacity in the power law above, so \$\$c_R(w) =
\kappa\\ w^{-\lambda}\$\$ for all \\w\\ less than `w_pp_cutoff` and zero
for larger sizes. Changing `kappa` therefore rescales the carrying
capacity. It has a second role in that the same expression also set the
initial resource abundance when the model was created: \$\$N_R(w) =
\kappa\\ w^{-\lambda}.\$\$ Unlike the carrying capacity, however, the
initial resource abundance is **not** updated when you subsequently
change `kappa` (or call
[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md)).

Assigning to `resource_params` only rebuilds the size-dependent resource
rate and capacity arrays from these scalars (leaving any arrays you have
set manually untouched). It does **not** balance the resource, i.e. it
does not adjust one of the rate or capacity to keep the resource at the
steady state where it replenishes at the rate at which it is consumed.
This mirrors the way the species parameters feed the species rates. If
you want to preserve the steady state after changing a resource scalar,
call
[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md)
with the appropriate argument (which balances by default).

## See also

[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md)
