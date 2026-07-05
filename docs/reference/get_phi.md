# Get values from feeding kernel function

This involves finding the feeding kernel function for each species,
using the pred_kernel_type parameter in the species_params data frame,
checking that it is valid and all its arguments are contained in the
species_params data frame, and then calling this function with the ppmr
vector.

## Usage

``` r
get_phi(species_params, ppmr)
```

## Arguments

- species_params:

  A species parameter data frame

- ppmr:

  Values of the predator/prey mass ratio at which to evaluate the
  predation kernel function

## Value

An array (species x ppmr) with the values of the predation kernel
function
