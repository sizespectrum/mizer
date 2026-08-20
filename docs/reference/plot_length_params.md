# The weight-length parameters to plot each row of plotting data with

A length axis needs an allometric weight-length relationship for every
line on the plot. The species take theirs from their species parameters
and the resource takes its from
[`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md),
where it defaults to the equivalent spherical diameter (see
[`resource_length_params()`](https://sizespectrum.org/mizer/reference/resource_length_params.md)).
Anything else — the "Total" row, for instance — has none, and is
reported as `NA` so that the caller can leave it out.

## Usage

``` r
plot_length_params(species, params)
```

## Arguments

- species:

  A vector of the species names in the plotting data.

- params:

  A MizerParams object providing the weight-length parameters.

## Value

A data frame with columns `a` and `b`, one row for each element of
`species`, holding `NA` where no relationship is known.
