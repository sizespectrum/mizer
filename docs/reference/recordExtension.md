# Record an extension and its version stamp on a mizer object

Writes an entry for `name` into the object's extension metadata,
converting it to the versioned list form. Existing entries (and their
version stamps) are preserved, keeping their position in the chain. A
genuinely new entry is prepended to the front of the list.

## Usage

``` r
recordExtension(params, name, version = NULL, requirement = NULL)
```

## Arguments

- params:

  A `MizerParams` object.

- name:

  The extension identifier (package/class name).

- version:

  Optional version string to stamp. If `NULL` (default) the existing
  stamp is preserved.

- requirement:

  Optional requirement string (e.g. `"sizespectrum/mizerMR"` or
  `"1.0.0"`).

## Value

The `params` object with updated extension metadata.

## Details

This is infrastructure for extension package authors. An extension
constructor normally calls `recordExtension()` and then
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md).
Ordinary model users do not need either call.

## See also

"Creating a mizer extension package": [Creating a mizer extension
package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.html)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`other_mort()`](https://sizespectrum.org/mizer/reference/other_mort.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
