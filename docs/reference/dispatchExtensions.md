# Filter an extension vector to those that participate in S3/S4 dispatch

An extension participates in dispatch if its requirement is
`NA_character_` (in-development) or if an S4 class with its name already
exists.

## Usage

``` r
dispatchExtensions(extensions)
```

## Arguments

- extensions:

  Named character vector of extensions.

## Value

A named character vector containing only the dispatch extensions,
preserving order.
