# Load (and optionally install) namespaces for all non-NA extensions

Load (and optionally install) namespaces for all non-NA extensions

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
