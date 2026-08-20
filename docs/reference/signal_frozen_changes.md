# Signal the changes to species parameters that cannot take effect

Goes through the rate arrays that can be frozen, see
[`frozen_rate_params()`](https://sizespectrum.org/mizer/reference/frozen_rate_params.md),
and raises a
[`signal_frozen()`](https://sizespectrum.org/mizer/reference/signal_frozen.md)
condition for each frozen array that one of the changed species
parameters feeds. This is what turns "the model no longer follows the
species parameters" into a warning the user sees at the moment they make
the change. It is one of the diagnostics that only
[`given_species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
gives, see there.

## Usage

``` r
signal_frozen_changes(params, changed)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object, holding the rate arrays as they are, that is, before the
  change is applied.

- changed:

  A character vector with the names of the species parameters that the
  user changed.

## Value

`NULL` invisibly. Called for its side effect of signalling.
