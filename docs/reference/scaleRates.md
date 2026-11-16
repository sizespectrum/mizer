# Rescale all rates in a mizer model

**\[experimental\]** Multiplies all rates in the model by a given
factor. Rescaling all rates by a factor \\f\\ is equivalent to rescaling
time by \\f\\: it speeds up (or slows down) all dynamics without
affecting the steady state of each species and resource. Extension
providers must declare how their dynamics transform under this change of
time scale.

## Usage

``` r
scaleRates(params, factor, ...)
```

## Arguments

- params:

  A MizerParams object

- factor:

  The positive factor by which all rates are multiplied.

- ...:

  Currently unused.

## Value

The MizerParams object with all rates rescaled by `factor`.

## Details

The following rates and their associated species parameters are
rescaled:

- Search volume (`search_vol` slot and `gamma` species parameter)

- Maximum intake rate (`intake_max` slot and `h` species parameter)

- Metabolic rate (`metab` slot and `ks`, `k` species parameters)

- External mortality (`mu_b` slot and `z0`, `z_ext`, `z0pre` species
  parameters)

- External encounter rate (`ext_encounter` slot and `E_ext` species
  parameter)

- External diffusion (`ext_diffusion` slot and `D_ext` species
  parameter)

- Catchability (`catchability` slot and `catchability` column in
  `gear_params`)

- Reproduction rate parameters (`R_max` and `constant_reproduction`
  multiply by the factor; `ricker_b` and `sheperd_b` divide by it)

- Resource growth rate (`rr_pp` slot and `r_pp` resource parameter)

Both the rate arrays stored in the MizerParams slots and the associated
species parameters in `species_params` and `given_species_params` are
rescaled, so that the parameters remain consistent with the rate arrays.
Changed scalar coefficients are recorded as given values, including
those previously calculated by mizer. Existing manually set arrays
retain their protection from recalculation.

Components and custom providers require declared transformation support.
Their hooks run after the core fields are scaled and before final
validation, with the original model available for reference. Component
states are left to their hooks; there is no automatic state rescaling.

## See also

[`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.md)
