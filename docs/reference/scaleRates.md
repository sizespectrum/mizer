# Rescale all rates in a mizer model

**\[experimental\]** Multiplies all rates in the model by a given
factor. Rescaling all rates by a factor \\f\\ is equivalent to rescaling
time by \\f\\: it speeds up (or slows down) all dynamics without
affecting the steady state of each species, provided the resource
spectrum is held at its steady-state value.

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

- Maximum reproduction rate (`R_max` species parameter)

- Resource growth rate (`rr_pp` slot)

Both the rate arrays stored in the MizerParams slots and the associated
species parameters in `species_params` and `given_species_params` are
rescaled, so that the parameters remain consistent with the rate arrays.

## See also

[`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.md)
