# Gaussian-mixture predation kernel

A predation kernel for which the log predator/prey mass ratio follows a
mixture of Gaussian distributions.

## Usage

``` r
gaussian_mixture_pred_kernel(ppmr, kernel_p, kernel_mean, kernel_sd)
```

## Arguments

- ppmr:

  A vector of predator/prey mass ratios.

- kernel_p:

  A numeric vector of relative component proportions.

- kernel_mean:

  A numeric vector of component means on the log predator/prey
  mass-ratio scale.

- kernel_sd:

  A numeric vector of positive component standard deviations.

## Value

A vector giving the value of the predation kernel at each of the
predator/prey mass ratios in the `ppmr` argument.

## Details

Writing the predator mass as \\w\\, the prey mass as \\w_p\\, and \\x =
\ln(w / w_p)\\, the feeding kernel is \$\$ \phi_i(w, w_p) = \sum_j
a\_{ij} \exp\left\[-\frac{(x - \mu\_{ij})^2}{2\sigma\_{ij}^2}\right\],
\qquad a\_{ij} = \frac{p\_{ij}/\sigma\_{ij}} {\sum_k
p\_{ik}/\sigma\_{ik}}. \$\$ for predator/prey mass ratios greater than
or equal to one, and zero for smaller ratios.

This is proportional to the Gaussian-mixture probability density with
mixing proportions \\p\_{ij}\\, means \\\mu\_{ij}\\, and standard
deviations \\\sigma\_{ij}\\. The scaling makes the sum of the component
peak heights equal to one. Consequently the kernel is at most one, and a
one-component mixture is identical to
[`lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/lognormal_pred_kernel.md)
with `beta = exp(kernel_mean)` and `sigma = kernel_sd`.

The three component parameters are vectors of equal length. When this
function is selected in a species parameter data frame, they should be
held in the list-columns `kernel_p`, `kernel_mean`, and `kernel_sd`. The
values in `kernel_p` must be non-negative with at least one positive
value, but they do not need to sum to one because they are normalised by
the function.

## See also

[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)

Other predation kernel:
[`box_pred_kernel()`](https://sizespectrum.org/mizer/reference/box_pred_kernel.md),
[`lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/lognormal_pred_kernel.md),
[`power_law_pred_kernel()`](https://sizespectrum.org/mizer/reference/power_law_pred_kernel.md),
[`truncated_lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/truncated_lognormal_pred_kernel.md)

## Examples

``` r
ppmr <- exp(seq(0, 12, length.out = 200))
phi <- gaussian_mixture_pred_kernel(
    ppmr,
    kernel_p = c(0.3, 0.7),
    kernel_mean = c(4, 8),
    kernel_sd = c(0.8, 1.5)
)
plot(ppmr, phi, type = "l", log = "x")
```
