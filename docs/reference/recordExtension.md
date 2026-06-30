# Record an extension and its version stamp on a mizer object

Writes an entry for `name` into the object's `@extensions` slot,
converting the slot to the versioned list form. Existing entries (and
their version stamps) are preserved. The requirement is taken from the
existing entry if present, otherwise from the registered extension
chain.

## Usage

``` r
recordExtension(params, name, version = NULL)
```

## Arguments

- params:

  A `MizerParams` object.

- name:

  The extension identifier (its S4 marker class name).

- version:

  Optional version string to stamp. If `NULL` (default) the existing
  stamp is preserved.

## Value

The `params` object with the updated `@extensions` slot.

## Details

Extension packages should call this instead of assigning to
`@extensions` directly. Pass `version` (typically
`packageVersion(name)`) only when the object has just been created or
upgraded to conform to that version; leave it `NULL` for ordinary
modifications so the existing stamp is preserved.

## See also

"Creating a mizer extension package":
[`vignette("creating-extension-packages", package = "mizer")`](https://sizespectrum.org/mizer/articles/creating-extension-packages.md)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`clearExtensionChain()`](https://sizespectrum.org/mizer/reference/clearExtensionChain.md),
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md),
[`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
