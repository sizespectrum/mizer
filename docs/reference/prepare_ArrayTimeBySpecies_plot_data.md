# The complete plotting data of a time-by-species array

The species selection and the background grouping are done in a single
pass, so that no species can be both selected under its own name and
appended again under the `"Background"` legend — which is what a
separate appending step used to do to every background species whenever
`species` was left at its default of all of them.

## Usage

``` r
prepare_ArrayTimeBySpecies_plot_data(
  x,
  species = NULL,
  tlim = c(NA, NA),
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  log_y = TRUE
)
```

## Arguments

- x:

  An `ArrayTimeBySpecies` object.

- species:

  Character vector of species to include, or `NULL` for all.

- tlim:

  Numeric vector of length two giving the time limits.

- ylim:

  Numeric vector of length two giving the value limits. Values outside
  them are dropped, as they cannot be seen anyway.

- total:

  Whether to append the total, which is the total over every species the
  array holds, whatever is drawn.

- background:

  Whether background species are included.

- log_y:

  Whether the values will be drawn on a logarithmic axis. Only then are
  non-positive values dropped: they have no place on a log axis, but on
  a linear one they are data like any other, and a quantity that can go
  negative — a rate of change, a difference between two models — would
  otherwise lose exactly the part of it that is interesting.

## Value

A data frame with `Year`, value, `Species` and `Legend` columns.
