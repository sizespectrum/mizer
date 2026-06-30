# Extract the requirement view of an extension chain

The `@extensions` slot may be stored either as a named character vector
of requirement strings (the legacy/unversioned form) or as a named list
whose entries are length-2 character vectors
`c(requirement = ..., version = ...)`. This helper returns the
requirement strings as a plain named character vector, which is the form
used for dispatch and suffix comparison.

## Usage

``` r
extensionRequirements(ext)
```

## Arguments

- ext:

  The contents of an `@extensions` slot (character vector or list).

## Value

A named character vector of requirement strings.
