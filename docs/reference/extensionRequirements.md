# Extract the requirement view of an extension chain

The `$extensions` element may be stored either as a named character
vector of requirement strings (the legacy/unversioned form) or as a
named list whose entries are length-2 character vectors
`c(requirement = ..., version = ...)`. This helper returns the
requirement strings as a plain named character vector.

## Usage

``` r
extensionRequirements(ext)
```

## Arguments

- ext:

  The contents of an `$extensions` element (character vector or list).

## Value

A named character vector of requirement strings.
