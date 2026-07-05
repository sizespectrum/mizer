# Set defaults for predation kernel parameters

If the predation kernel type has not been specified for a species, then
it is set to "lognormal" and the default values are set for the
parameters `beta` and `sigma`.

## Usage

``` r
default_pred_kernel_params(object)
```

## Arguments

- object:

  Either a MizerParams object or a species parameter data frame

## Value

The `object` with updated columns in the species params data frame.
