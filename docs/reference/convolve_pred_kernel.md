# Sum over the predation kernel

**\[experimental\]** Evaluates the sums over the predation kernel that
underlie the encounter rate, the predation rate and the predation
diffusion, using the Fourier-transformed kernels that
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
stores in the MizerParams object. This is how
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md),
[`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md),
[`mizerDiffusion()`](https://sizespectrum.org/mizer/reference/mizerDiffusion.md)
and [`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md)
evaluate these sums. An extension package that needs such a sum for a
density of its own should call this function rather than multiplying the
Fourier-transformed kernels itself, so that it does not depend on how
mizer stores them.

## Usage

``` r
convolve_pred_kernel(
  params,
  x,
  kernel = c("encounter", "predation", "diffusion")
)
```

## Arguments

- params:

  A MizerParams object whose predation kernel depends only on the
  predator/prey mass ratio, i.e., one without a custom predation kernel
  set with
  [`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md).

- x:

  For `kernel = "encounter"` or `"diffusion"`, a matrix (species x size
  in `w_full(params)`) of prey densities, or a vector with one entry per
  size in `w_full(params)` that is used for every species. For
  `kernel = "predation"`, a matrix (species x size in `w(params)`) of
  predator densities, or a vector with one entry per size in
  `w(params)`.

- kernel:

  Which sum to evaluate: "encounter" (the default), "predation" or
  "diffusion".

## Value

A matrix (species x predator size in `w(params)`) for
`kernel = "encounter"` or `"diffusion"`, and a matrix (species x prey
size in `w_full(params)`) for `kernel = "predation"`.

## Details

With `kernel = "encounter"` the result is \$\$\sum\_{k} \phi_i(w_j,
w_k)\\ x\_{ik}\$\$ for each species \\i\\ and each predator size \\w_j\\
in `w(params)`, where the sum runs over the prey sizes \\w_k\\ in
`w_full(params)`. With `kernel = "diffusion"` the same sum is evaluated
with the kernel used for the predation diffusion, which differs from the
encounter kernel only when bin averaging is switched on (see
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)).
With `kernel = "predation"` the result is \$\$\sum\_{j} \phi_i(w_j,
w_k)\\ x\_{ij}\$\$ for each species \\i\\ and each prey size \\w_k\\ in
`w_full(params)`, where the sum runs over the predator sizes \\w_j\\ in
`w(params)`.

The entries of `x` are the summands without the kernel, so they include
the size-bin widths: for the encounter rate, for example, `x` is the
prey biomass density times `dw_full(params)`. Pair them with the plain
point weights, not with bin-averaged ones, even when bin averaging is
switched on; see
[`encounter_kernel()`](https://sizespectrum.org/mizer/reference/encounter_kernel.md)
for why.

The sums are evaluated as convolutions by fast Fourier transform. For a
model of defaults edition 3 or later (see
[`defaults_edition()`](https://sizespectrum.org/mizer/reference/defaults_edition.md))
they only run over sizes on the grid: the part of a kernel that extends
below the smallest size in `w_full(params)` contributes nothing. This is
achieved by evaluating the convolutions on a zero-padded grid, so that
the part of a kernel that falls off one end of the grid does not wrap
around to pick up the densities at the other end. A model of an earlier
edition keeps the circular convolutions of earlier versions of mizer, so
that its results do not change. In such a model, the part of the kernel
of a small predator that extends below the smallest size in
`w_full(params)` wraps around to the largest sizes, so that the smallest
predators are credited with eating the largest prey.

## See also

[`encounter_kernel()`](https://sizespectrum.org/mizer/reference/encounter_kernel.md),
[`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)

## Examples

``` r
params <- NS_params
# The rate at which each predator encounters resource biomass, per unit
# search volume
prey <- initialNResource(params) * w_full(params) * dw_full(params)
avail <- convolve_pred_kernel(params, prey)
# Add the interaction with the resource and the search volume to get the
# resource part of the encounter rate
encounter_resource <- avail * species_params(params)$interaction_resource *
    search_vol(params)
```
