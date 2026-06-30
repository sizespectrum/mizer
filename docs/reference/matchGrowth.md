# Adjust model to produce observed growth

**\[experimental\]** Scales the search volume, the maximum consumption
rate, the metabolic rate and the external encounter rate all by the same
factor in order to achieve a growth rate that allows individuals to
reach their maturity size by their maturity age while keeping the
feeding level and the critical feeding level unchanged. Then
recalculates the size spectra using
[`steadySingleSpecies()`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md).

## Usage

``` r
matchGrowth(params, species = NULL, keep = c("egg", "biomass", "number"), ...)
```

## Arguments

- params:

  A MizerParams object

- species:

  The species to be affected. Optional. By default all species for which
  growth information is available will be affected. A vector of species
  names, or a numeric vector with the species indices, or a logical
  vector indicating for each species whether it is to be affected (TRUE)
  or not.

- keep:

  A string determining which quantity is to be kept constant. The
  choices are "egg" which keeps the egg density constant, "biomass"
  which keeps the total biomass of the species constant and "number"
  which keeps the total number of individuals constant.

- ...:

  Additional arguments passed to the method.

## Value

A modified MizerParams object with rescaled search volume, maximum
consumption rate and metabolic rate and rescaled species parameters
`gamma`,`h`, `ks` and `k`.

## Details

Maturity size and age are taken from the `w_mat` and `age_mat` columns
in the species_params data frame. If `age_mat` is missing, mizer
calculates it from the von Bertalanffy growth curve parameters using
[`age_mat_vB()`](https://sizespectrum.org/mizer/reference/age_mat_vB.md).
If those are not available either for a species, the growth rate for
that species will not be changed.

## Examples

``` r
# Rescale rates so all species reach maturity by their maturity age.
# The search volume gamma is adjusted to achieve the correct growth rate.
species_params(NS_params)["Cod", "gamma"]
#> [1] 1.599016e-10
params <- matchGrowth(NS_params)
species_params(params)["Cod", "gamma"]
#> [1] 2.351462e-10
age_mat(params)["Cod"]
#>      Cod 
#> 1.954054 
```
