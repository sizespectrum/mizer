# Extract the version stamps of an extension chain

Returns the version of the extension package that last upgraded the
object for each extension, or `NA_character_` where no stamp is recorded
(including the legacy character-vector form, which carries no versions).

## Usage

``` r
extensionVersions(ext)
```

## Arguments

- ext:

  The contents of an `@extensions` slot (character vector or list).

## Value

A named character vector of version strings (`NA` where unknown).
