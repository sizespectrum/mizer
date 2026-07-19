# Set predation kernel

You will usually not need to call this function directly. Instead change
the relevant species parameters (`pred_kernel_type`, and, depending on
its value, `beta`/`sigma` or `ppmr_min`/`ppmr_max`) with
`given_species_params(params) <-` and let mizer recalculate the
predation kernel for you. Call `setPredKernel()` directly only if you
want to supply the full kernel array yourself. See
[`vignette("cheatsheet-changing-parameters")`](https://sizespectrum.org/mizer/articles/cheatsheet-changing-parameters.md)
for a full explanation of when to reach for which level of the model.

## Usage

``` r
setPredKernel(params, pred_kernel = NULL, reset = FALSE, ...)

getPredKernel(params)

pred_kernel(params)

pred_kernel(params) <- value
```

## Arguments

- params:

  A MizerParams object

- pred_kernel:

  Optional. An array (species x predator size x prey size) that holds
  the predation coefficient of each predator at size on each prey size.
  If not supplied, a default is set as described in section "Setting
  predation kernel".

- reset:

  If set to TRUE, then the predation kernel will be reset to the value
  calculated from the species parameters, even if it was previously
  overwritten with a custom value. If set to FALSE (default) then a
  recalculation from the species parameters will take place only if no
  custom value has been set.

- ...:

  Unused

- value:

  pred_kernel

## Value

`setPredKernel()`: A MizerParams object with updated predation kernel.

`getPredKernel()` or equivalently `pred_kernel()`: An array (predator
species x predator_size x prey_size)

## Details

The predation kernel determines the distribution of prey sizes that a
predator feeds on. It is used in
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
when calculating the rate at which food is encountered and in
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md)
when calculating the rate at which a prey is predated upon. The
predation kernel can be a function of the predator/prey size ratio or it
can be a function of the predator size and the prey size separately.
Both types can be set up with this function.

## Setting predation kernel

**Kernel dependent on predator to prey size ratio**

If the `pred_kernel` argument is not supplied, then this function sets a
predation kernel that depends only on the ratio of predator mass to prey
mass, not on the two masses independently. The shape of that kernel is
then determined by the `pred_kernel_type` column in species_params.

The default for `pred_kernel_type` is "lognormal". This will call the
function
[`lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/lognormal_pred_kernel.md)
to calculate the predation kernel. An alternative pred_kernel type is
"box", implemented by the function
[`box_pred_kernel()`](https://sizespectrum.org/mizer/reference/box_pred_kernel.md),
and "power_law", implemented by the function
[`power_law_pred_kernel()`](https://sizespectrum.org/mizer/reference/power_law_pred_kernel.md).
These functions require certain species parameters in the species_params
data frame. For the lognormal kernel these are `beta` and `sigma`, for
the box kernel they are `ppmr_min` and `ppmr_max`. They are explained in
the help pages for the kernel functions. Except for `beta` and `sigma`,
no defaults are set for these parameters. If they are missing from the
species_params data frame then mizer will issue an error message.

You can use any other string for `pred_kernel_type`. If for example you
choose "my" then you need to define a function `my_pred_kernel` that you
can model on the existing functions like
[`lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/lognormal_pred_kernel.md).

When using a kernel that depends on the predator/prey size ratio only,
mizer does not need to store the entire three dimensional array in the
MizerParams object. Such an array can be very big when there is a large
number of size bins. Instead, mizer only needs to store two
two-dimensional arrays that hold Fourier transforms of the feeding
kernel function that allow the encounter rate and the predation rate to
be calculated very efficiently. However, if you need the full
three-dimensional array you can calculate it with the `getPredKernel()`
function.

**Kernel dependent on both predator and prey size**

If you want to work with a feeding kernel that depends on predator mass
and prey mass independently, you can specify the full feeding kernel as
a three-dimensional array (predator species x predator size x prey
size).

You should use this option only if a kernel dependent only on the
predator/prey mass ratio is not appropriate. Using a kernel dependent on
predator/prey mass ratio only allows mizer to use fast Fourier transform
methods to significantly reduce the running time of simulations.

The order of the predator species in `pred_kernel` should be the same as
the order in the species params dataframe in the `params` object. If you
supply a named array then the function will check the order and warn if
it is different.

## Higher-order quadrature

When the predation kernel depends only on the predator/prey mass ratio,
the encounter and predation rates are evaluated as convolutions using
the fast Fourier transform. By default mizer point-samples the kernel at
the grid nodes, which is a first-order (rectangle-rule) quadrature. When
the `bin_average` entry of the `second_order_w` slot is `TRUE`, the
kernel is instead integrated over each logarithmic size bin when
building the Fourier-transformed kernels. This finite-volume consistent
quadrature lifts the encounter and predation rates towards second order
at no extra runtime cost, because the integration is performed once here
and the rate functions themselves are unchanged. The predation kernel is
additionally averaged over the prey bin (a trapezoid fold), so that the
predation rate returned by
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md)
is the prey-bin average that the predation-mortality sink needs to be
second order. The default remains the first-order scheme so that
existing models are unaffected. Enable it with
`second_order_w(params) <- TRUE` (see
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)),
which also turns on the other bin-averaged rate quadratures so the whole
model stays consistent. See the
[`vignette("fft")`](https://sizespectrum.org/mizer/articles/fft.md) for
the mathematical details.

## See also

Other functions for setting parameters:
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
[`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md),
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md),
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)

## Examples

``` r
## Set up a MizerParams object
params <-  NS_params

## If you change predation kernel parameters after setting up a model,
#  this will be used to recalculate the kernel
species_params(params)["Cod", "beta"] <- 200

## You can change to a different predation kernel type
species_params(params)$ppmr_max <- 4000
species_params(params)$ppmr_min <- 200
species_params(params)$pred_kernel_type <- "box"
plot(w_full(params), getPredKernel(params)["Cod", 100, ], type="l", log="x")


## If you need a kernel that depends also on prey size you need to define
# it yourself.
pred_kernel <- getPredKernel(params)
pred_kernel["Herring", , ] <- sweep(pred_kernel["Herring", , ], 2,
                                    params@w_full, "*")
params<- setPredKernel(params, pred_kernel = pred_kernel)
```
