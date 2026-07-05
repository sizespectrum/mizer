# Register a single mizer extension for this R session

Prepends one extension to the front of the active extension chain,
giving it the highest dispatch priority. Designed to be called from a
package's `.onLoad` hook so that the chain grows naturally in load
order: the last package loaded ends up outermost.

## Usage

``` r
registerExtension(name, requirement = NA_character_, install = FALSE)
```

## Arguments

- name:

  A syntactically valid R name identifying the extension (e.g.
  `"mizerExtA"`). This name is used as the S4 marker class name.

- requirement:

  A version string, installation specification, or `NA_character_` (the
  default). `NA_character_` marks an in-development extension whose S4
  marker class mizer creates automatically. A version string such as
  `"1.2.0"` records the minimum required package version.

- install:

  Logical. If `TRUE`, attempt to install a missing extension package.

## Value

The updated extension chain, invisibly.

## Details

The call is idempotent: if the extension is already registered at any
position in the chain, the function returns silently without modifying
the chain. This makes it safe to call from
[`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html),
which re-executes `.onLoad`.

## See also

[`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md)
for registering an explicit full chain. "Using mizer extension
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
[`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
