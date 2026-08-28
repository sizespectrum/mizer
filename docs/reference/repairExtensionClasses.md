# Rebuild the marker classes of an extension chain if any went missing

[`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
removes the S4 classes held by the namespace it reloads, which can take
a marker class out of the chain. Recreating only the missing class is
not enough: R prunes a removed superclass from the `contains` list of
its subclasses, so a marker class that used to sit outside the missing
one is left parented directly on `MizerParams` and would make
[`defineOrCheckClass()`](https://sizespectrum.org/mizer/reference/defineOrCheckClass.md)
stop. The whole dynamic chain is therefore removed and rebuilt,
outermost first so that no class is removed while a subclass of it still
exists.

## Usage

``` r
repairExtensionClasses(extensions)
```

## Arguments

- extensions:

  Named character vector of extensions (full chain or dispatch subset).

## Value

Invisibly, `TRUE` if a repair was carried out, `FALSE` otherwise.

## Details

The chain is inspected first and left completely untouched when it is
intact, which is the usual case on the repeated-registration path.
