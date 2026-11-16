# Determine reproduction rate needed for initial egg abundance

Determine reproduction rate needed for initial egg abundance

## Usage

``` r
getRequiredRDD(params, ..., rates = NULL)
```

## Arguments

- params:

  A MizerParams object

- ...:

  Arguments passed to methods.

- rates:

  Optional already evaluated rate list for this model's initial state
  and effort, used to share an effective evaluation with a caller. Must
  contain growth, total mortality and diffusion. Usually omitted.

## Value

A vector of reproduction rates for all species
