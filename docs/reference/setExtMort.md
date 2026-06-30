# Set external mortality rate

Set external mortality rate

## Usage

``` r
setExtMort(
  params,
  ext_mort = NULL,
  z0pre = 0.6,
  z0exp = params@resource_params$n - 1,
  reset = FALSE,
  z0 = deprecated(),
  ...
)

getExtMort(params)

ext_mort(params)

ext_mort(params) <- value
```

## Arguments

- params:

  MizerParams

- ext_mort:

  Optional. An array (species x size) holding the external mortality
  rate. If not supplied, a default is set as described in the section
  "Setting external mortality rate".

- z0pre:

  If `z0`, the mortality from other sources, is not a column in the
  species data frame, it is calculated as z0pre \* w_inf ^ z0exp.
  Default value is 0.6.

- z0exp:

  If `z0`, the mortality from other sources, is not a column in the
  species data frame, it is calculated as `z0pre * w_inf ^ z0exp`.
  Default value is `n-1`.

- reset:

  If set to TRUE, then the external mortality rate will be reset to the
  value calculated from the species parameters, even if it was
  previously overwritten with a custom value. If set to FALSE (default)
  then a recalculation from the species parameters will take place only
  if no custom value has been set.

- z0:

  **\[deprecated\]** Use `ext_mort` instead. Not to be confused with the
  species_parameter `z0`.

- ...:

  Unused

- value:

  ext_mort

## Value

`setExtMort()`: A MizerParams object with updated external mortality
rate.

`getExtMort()` or equivalently `ext_mort()`: An `ArraySpeciesBySize`
object (species x size) with the external mortality.

## Setting external mortality rate

The external mortality is all the mortality that is not due to fishing
or predation by predators included in the model. The external mortality
could be due to predation by predators that are not explicitly included
in the model (e.g. mammals or seabirds) or due to other causes like
illness. It is a rate with units 1/year.

The `ext_mort` argument allows you to specify an external mortality rate
that depends on species and body size. You can see an example of this in
the Examples section of the help page for `setExtMort()`.

If the `ext_mort` argument is not supplied, then the external mortality
is taken from the species parameters as \$\$\mu\_{ext.i}(w) = z\_{0.i} +
z\_{ext.i} w^{d_i}.\$\$ The value of the constant \\z_0\\ for each
species is taken from the `z0` column of the species parameter data
frame, if that column exists. Otherwise it is calculated as \$\$z\_{0.i}
= {\tt z0pre}\_i\\ w\_{inf}^{\tt z0exp}.\$\$ Missing values of `z_ext`
are set to 0 and missing values of `d` are set to `n - 1`.

By default the power law is evaluated at the left bin edges \\w_j\\
(point sampling). If the `bin_average` entry of the `second_order_w`
slot is `TRUE` (see
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)),
then the \\z\_{ext} w^d\\ term is instead replaced by its exact average
over each bin \\\[w_j, w\_{j+1}\]\\, \$\$\frac{z\_{ext}}{\Delta
w_j}\int\_{w_j}^{w\_{j+1}} w^d\\ dw = z\_{ext}\\\frac{w\_{j+1}^{d+1} -
w_j^{d+1}}{(d+1)\\\Delta w_j},\$\$ (with the limiting form
\\z\_{ext}\ln(w\_{j+1}/w_j)/\Delta w_j\\ when \\d = -1\\). This is the
consistent choice in the finite-volume scheme, where the external
mortality multiplies the bin-averaged abundance. The bin-averaging is
applied only to the auto-calculated power-law default; a user-supplied
`ext_mort` array is left untouched.

## See also

Other functions for setting parameters:
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
[`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md),
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md),
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)

## Examples

``` r
params <- newMultispeciesParams(NS_species_params)
#> Because you have n != p, the default value for `h` is not very good.
#> Because the age at maturity is not known, I need to fall back to using
#> von Bertalanffy parameters, where available, and this is not reliable.
#> No ks column so calculating from critical feeding level.
#> Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.
#> Using f0, h, lambda, kappa and the predation kernel to calculate gamma.

#### Setting allometric death rate #######################

# Set coefficient for each species. Here we choose 0.1 for each species
z0pre <- rep(0.1, nrow(species_params(params)))

# Multiply by power of size with exponent, here chosen to be -1/4
# The outer() function makes it an array species x size
allo_mort <- outer(z0pre, w(params)^(-1/4))

# Change the external mortality rate in the params object
ext_mort(params) <- allo_mort
```
