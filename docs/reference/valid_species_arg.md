# Helper function to assure validity of species argument

If the species argument contains invalid species, then these are ignored
but a warning is issued.

## Usage

``` r
valid_species_arg(
  object,
  species = NULL,
  return.logical = FALSE,
  error_on_empty = FALSE
)
```

## Arguments

- object:

  A MizerSim or MizerParams object from which the species should be
  selected.

- species:

  The species to be selected. Optional. By default all target species
  are selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- return.logical:

  Whether the return value should be a logical vector. Default FALSE.

- error_on_empty:

  Whether to throw an error if there are zero valid species. Default
  FALSE.

## Value

A vector of species names, in the same order as specified in the
'species' argument. If 'return.logical = TRUE' then a logical vector is
returned instead, with length equal to the number of species, with TRUE
entry for each selected species.
