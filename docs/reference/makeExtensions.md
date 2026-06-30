# Build a versioned extension list from requirements and versions

Build a versioned extension list from requirements and versions

## Usage

``` r
makeExtensions(requirements, versions = character())
```

## Arguments

- requirements:

  A named character vector of requirement strings.

- versions:

  A named character vector of version stamps. Names not present default
  to `NA_character_`.

## Value

A named list whose entries are `c(requirement = ..., version = ...)`, or
an empty character vector when `requirements` is empty.
