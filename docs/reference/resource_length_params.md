# The weight-length parameters of the resource

Reads `a` and `b` from
[`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md),
falling back to
[resource_length_defaults](https://sizespectrum.org/mizer/reference/resource_length_defaults.md)
for a model that does not set them — which is every model built before
these parameters existed.

## Usage

``` r
resource_length_params(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

A list with entries `a` and `b`.
