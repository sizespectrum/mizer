# Define S4 marker classes for a set of dispatch extensions

Creates a linear inheritance chain of S4 classes: the outermost
extension extends the next, which extends the next, down to the base
`MizerParams` / `MizerSim` class. Existing classes are checked for
compatibility instead of being redefined.

## Usage

``` r
defineExtensionClasses(extensions)
```

## Arguments

- extensions:

  Named character vector of extensions (full chain or dispatch subset).
  Non-dispatch entries are silently ignored.

## Value

Invisibly, the named character vector of dispatch extensions.
