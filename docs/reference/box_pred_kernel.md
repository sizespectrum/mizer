# Box predation kernel

A predation kernel where the predator/prey mass ratio is uniformly
distributed on an interval.

## Usage

``` r
box_pred_kernel(ppmr, ppmr_min, ppmr_max)
```

## Arguments

- ppmr:

  A vector of predator/prey size ratios

- ppmr_min:

  Minimum predator/prey mass ratio

- ppmr_max:

  Maximum predator/prey mass ratio

## Value

A vector giving the value of the predation kernel at each of the
predator/prey mass ratios in the `ppmr` argument.

## Details

Writing the predator mass as \\w\\ and the prey mass as \\w_p\\, the
feeding kernel is 1 if \\w/w_p\\ is between `ppmr_min` and `ppmr_max`
inclusive and zero otherwise. `ppmr_min` must be strictly smaller than
`ppmr_max`. The parameters need to be given in the species parameter
dataframe in the columns `ppmr_min` and `ppmr_max`.

## See also

[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)

Other predation kernel:
[`lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/lognormal_pred_kernel.md),
[`power_law_pred_kernel()`](https://sizespectrum.org/mizer/reference/power_law_pred_kernel.md),
[`truncated_lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/truncated_lognormal_pred_kernel.md)

## Examples

``` r
params <- NS_params
# Set all required paramters before changing kernel type
species_params(params)$ppmr_max <- 4000
species_params(params)$ppmr_min <- 200
species_params(params)$pred_kernel_type <- "box"
plot(w_full(params), getPredKernel(params)["Cod", 10, ], type="l", log="x")
```
