# Load (and optionally install) namespaces for all non-NA extensions

For each extension whose requirement is not `NA_character_`, checks that
the package is installed and up-to-date, installs or upgrades via
[`pak::pkg_install()`](https://pak.r-lib.org/reference/pkg_install.html)
if `install = TRUE`, then calls
[`loadNamespace()`](https://rdrr.io/r/base/ns-load.html).

## Usage

``` r
ensureExtensionNamespaces(extensions, install = FALSE)
```

## Arguments

- extensions:

  Named character vector of extensions.

- install:

  Logical. If `TRUE`, install or upgrade missing/outdated packages via
  [`pak::pkg_install()`](https://pak.r-lib.org/reference/pkg_install.html).

## Value

Invisibly `TRUE`.
