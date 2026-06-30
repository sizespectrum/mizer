# Register mizer extensions for this R session

Registers an explicit full extension chain for the current R session.
The order of `extensions` is the S3 dispatch order, from outermost to
innermost extension. For example
`c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")` dispatches to `mizerExtB`
methods first, then `mizerExtA` methods, then base mizer methods.

## Usage

``` r
registerExtensions(extensions, install = FALSE)
```

## Arguments

- extensions:

  A named character vector. Names are extension identifiers. Values are
  version strings, installation specifications, or `NA_character_`.
  Installed extensions only participate in S3 dispatch if they provide
  an S4 marker class with the same name. `NA_character_` entries are
  treated as in-development dispatch extensions and mizer creates their
  marker classes automatically.

- install:

  Logical. If `TRUE`, missing or outdated extension packages are
  installed via
  [`pak::pkg_install()`](https://pak.r-lib.org/reference/pkg_install.html).
  Version strings install from CRAN; other requirement strings (e.g.
  `"user/repo@v1.2.0"`) are passed directly to pak and may refer to
  GitHub, local paths, or any other pak-supported source.

## Value

The active maximal extension chain, invisibly.

## Details

A session can handle objects whose extension chain is a suffix of the
registered maximal chain. For example, after registering
`c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")`, objects using only
`c(mizerExtA = "0.4.1")` are also valid.

For extension packages that register themselves incrementally from
`.onLoad`, use
[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md)
instead.

## See also

[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md)
for the incremental per-package variant. "Using mizer extension
packages":
[`vignette("using-extension-packages", package = "mizer")`](https://sizespectrum.org/mizer/articles/using-extension-packages.md).
"Creating a mizer extension package":
[`vignette("creating-extension-packages", package = "mizer")`](https://sizespectrum.org/mizer/articles/creating-extension-packages.md)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`clearExtensionChain()`](https://sizespectrum.org/mizer/reference/clearExtensionChain.md),
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
