# Signal the changes that are ignored because another parameter was given

Some species parameters are only used to calculate a default for another
one: `f0` for `gamma`, `fc` for `ks`, `age_mat` for `h`, and `k_vb` for
`h` or `age_mat`. Once the other one has been given, the model no longer
consults them, so changing them has no effect. This raises a warning
about that. It is one of the diagnostics that only
[`given_species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
gives, see there.

## Usage

``` r
signal_ignored_changes(given, changed)
```

## Arguments

- given:

  The given species parameters, as they are before the change.

- changed:

  A named list with one logical vector per changed column, saying which
  species were given a value, as built by
  [`given_species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md).

## Value

`NULL` invisibly. Called for its side effect of signalling.

## Details

Only a value that is there can be ignored, so this is asked about the
species that were *given* a value, not about every species whose value
changed: clearing a value to `NA` is a change, but not one this has
anything to say about.
