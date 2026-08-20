# Signal a gear parameter changed through the given species parameters

Mizer looks for the gear parameters in the gear parameter table, which
is read only when the model is built, so changing one of them through
the species parameters does not reach the model. This is one of the
diagnostics that only
[`given_species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
gives;
[`species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
stays quiet, see there.

## Usage

``` r
signal_gear_params_changes(changed)
```

## Arguments

- changed:

  A named list with one entry per changed column, or a character vector
  of the changed column names.

## Value

`NULL` invisibly. Called for its side effect of signalling.

## Details

`yield_observed` is the exception. It belongs in
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
which gives the observed yield per gear and species, and that is where
it should be set, which is what this reports. But it feeds no rate, so a
value in the species parameters is not lost:
[`get_yield_observed()`](https://sizespectrum.org/mizer/reference/get_yield_observed.md)
falls back to it for any species that has no observation among the gear
parameters.
