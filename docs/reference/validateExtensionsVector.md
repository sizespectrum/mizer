# Validate and normalise an extensions named character vector

Checks that `extensions` is a named character vector with unique,
syntactically valid names, and normalises `NULL` to
[`character()`](https://rdrr.io/r/base/character.html).

## Usage

``` r
validateExtensionsVector(extensions)
```

## Arguments

- extensions:

  A named character vector, or `NULL`.

## Value

A validated named character vector (possibly length-zero).
