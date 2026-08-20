# Calibrate the model scale to match a total observation

Internal implementation shared by
[`calibrateBiomass()`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md)
and
[`calibrateNumber()`](https://sizespectrum.org/mizer/reference/calibrateNumber.md).
Rescales the model with
[`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.md)
so that the total over all observed species of the modelled quantity
agrees with the total of the observations. Species with no observation
are left out of both totals.

## Usage

``` r
calibrate_to(params, to = c("biomass", "number"))
```

## Arguments

- params:

  A MizerParams object.

- to:

  The type of observation, either "biomass" or "number".

## Value

A MizerParams object. If no non-missing observations are provided, the
original object is returned unchanged.
