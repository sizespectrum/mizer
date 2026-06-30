# Upgrade the core slots of a MizerParams object

This is the `MizerParams` method of the
[`upgrade()`](https://rdrr.io/r/utils/upgrade.html) generic and performs
the *core mizer* upgrade only. It is called from the orchestrator
[`runExtensionUpgrades()`](https://sizespectrum.org/mizer/reference/runExtensionUpgrades.md),
which also invokes any registered extension upgrade methods. You should
never need to call it directly; use
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
(or
[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md))
instead.

## Usage

``` r
# S3 method for class 'MizerParams'
upgrade(object, ...)
```

## Arguments

- object:

  An old MizerParams object to be upgraded

- ...:

  Unused.

## Value

The upgraded MizerParams object

## See also

[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
