# Is a marker class one that mizer created dynamically?

Marker classes that mizer creates in
[`defineExtensionClasses()`](https://sizespectrum.org/mizer/reference/defineExtensionClasses.md)
belong to `.GlobalEnv` in the S4 class registry, although their metadata
bindings live in mizer's persistent class environment (see
[`extensionClassEnvironment()`](https://sizespectrum.org/mizer/reference/extensionClassEnvironment.md)).
This means mizer is free to remove and rebuild them. A class of the same
name that an extension package defines statically belongs to that
package's namespace and must never be touched.

## Usage

``` r
isDynamicMarkerClass(class)
```

## Arguments

- class:

  Character string — the S4 class name to test.

## Value

`TRUE` if `class` exists and carries the `.GlobalEnv` package identity
that marks it as one of mizer's dynamic classes.
