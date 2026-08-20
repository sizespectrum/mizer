# The modelled counterpart of an observation

Internal helper for the calibration and matching functions. Integrates
the initial abundance of each species over the sizes that the
observation covers and returns the total biomass in grams for
`to = "biomass"` or the total number of individuals for `to = "number"`.
The `<to>_cutoff` species parameter sets the smallest size counted for
each species; where it is missing the whole size range of the species is
counted.

## Usage

``` r
model_observation(params, to = c("biomass", "number"))
```

## Arguments

- params:

  A MizerParams object.

- to:

  The type of observation, either "biomass" or "number".

## Value

A named vector with one value for each species.

## Details

The integral is done by
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md),
so it follows the model's quadrature scheme and lets the bin straddling
the cutoff contribute only the part of it that lies above the cutoff.
