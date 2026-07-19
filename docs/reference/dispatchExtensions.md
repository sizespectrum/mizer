# Filter an extension vector to those that participate in S3/S4 dispatch

An extension participates in dispatch if its requirement is
`NA_character_` (in-development), if an S4 class with its name already
exists, or if its loaded package registers S3 dispatch methods for its
marker class (see
[`providesDispatchMethods()`](https://sizespectrum.org/mizer/reference/providesDispatchMethods.md)).
The last case lets an installed extension package participate without
defining its marker class statically; mizer creates the class
dynamically in
[`defineExtensionClasses()`](https://sizespectrum.org/mizer/reference/defineExtensionClasses.md).

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
