# Coerce a mizer object to its registered extension class

Coerces a `MizerParams` or `MizerSim` object to the S4 marker class
corresponding to the object's own extension chain. For `MizerSim`, the
extension chain is read from `sim@params@extensions`.

## Usage

``` r
coerceToExtensionClass(object, extensions = objectExtensions(object))
```

## Arguments

- object:

  A `MizerParams` or `MizerSim` object.

- extensions:

  Optional extension chain. Defaults to the chain stored in `object`, or
  in `object@params` for `MizerSim`.

## Value

The same object coerced to the appropriate marker class, or to the base
class for an empty extension chain.

## See also

"Creating a mizer extension package":
[`vignette("creating-extension-packages", package = "mizer")`](https://sizespectrum.org/mizer/articles/creating-extension-packages.md)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`clearExtensionChain()`](https://sizespectrum.org/mizer/reference/clearExtensionChain.md),
[`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md),
[`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
