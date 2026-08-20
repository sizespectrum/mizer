# The predation kernel as used by the encounter quadrature

**\[experimental\]** Returns the kernel array \\\Phi_i(w_k, w_p)\\ for
which \$\$E_i(w_k) = \gamma_i(w_k) \sum_p \Phi_i(w_k, w_p)
N^{eff}\_i(w_p) w_p \Delta w_p\$\$ reproduces exactly the available
energy computed by
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md),
where \\N^{eff}\\ is the interaction-weighted prey density. It is the
kernel that any summary function must use if its result is to be
consistent with
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md).

## Usage

``` r
encounter_kernel(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

An array (predator species x predator size x prey size).

## Details

On the default first-order path this is just the point-sampled kernel
returned by
[`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md).
When second-order bin-averaging is switched on (see
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md))
the two differ:
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
then builds the Fourier-transformed kernel from the kernel *integrated
over the prey bin*, divided by \\\beta - 1\\ so that the plain point
weight \\w_p \Delta w_p\\ carried by the prey vector is cancelled. Those
bin-integrated weights are recovered here from `params@ft_pred_kernel_e`
by an inverse Fourier transform, which costs one FFT and keeps this
helper automatically in step with whatever quadrature
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
used.

Pair it with the **plain point prey weight**
`params@w_full * params@dw_full`. That weight is a normalisation which
the kernel construction is built to cancel, not a first-order quadrature
weight, so it must not be passed through
[`bin_average_weight()`](https://sizespectrum.org/mizer/reference/bin_average_weight.md):
doing so applies the prey-bin integral twice. A summary function that
instead pairs the point-sampled
[`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
with a bin-averaged prey weight double-counts that quadrature; that was
the bug behind issue \#474.

## See also

[`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
for the point-sampled kernel used for plotting and for supplying a
custom kernel,
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md),
[`bin_average_weight()`](https://sizespectrum.org/mizer/reference/bin_average_weight.md)
