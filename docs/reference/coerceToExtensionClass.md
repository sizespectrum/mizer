# Coerce a mizer object to its extension class

Coerces a `MizerParams` or `MizerSim` object to the S3 class hierarchy
corresponding to the object's own extension list. For `MizerSim`, the
extension names are read from `sim$params$extensions` (or
`sim@params@extensions`).

## Usage

``` r
coerceToExtensionClass(object, extensions = objectExtensions(object))
```

## Arguments

- object:

  A `MizerParams` or `MizerSim` object.

- extensions:

  Optional extension list or vector. Defaults to the extensions stored
  in `object`, or in `object$params` for `MizerSim`.

## Value

The same object coerced to the appropriate S3 class vector, or to the
base class for an empty extension list.

## See also

"Creating a mizer extension package": [Creating a mizer extension
package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.html)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`other_mort()`](https://sizespectrum.org/mizer/reference/other_mort.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
