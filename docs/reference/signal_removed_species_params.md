# Signal that species parameter columns have been removed

A column that is absent from the table assigned to
[`species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
or
[`given_species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md)
is one the user no longer supplies. Mizer calculates afresh those it
knows how to calculate, and the rest leave the model altogether. This
reports the ones that leave the given species parameters, whether or not
mizer puts a calculated value back.

## Usage

``` r
signal_removed_species_params(withdrawn)
```

## Arguments

- withdrawn:

  A character vector of the withdrawn column names. Nothing is signalled
  when it is empty.

## Value

`NULL` invisibly. Called for its side effect of signalling.

## Details

The report is made at level 3, so it is visible in ordinary use and
silent from `info_level` 1 downwards: removing a column is what the user
asked for, not something that went differently from how they asked. It
carries the class `info_about_removed` for code that wants to catch it
in particular.
