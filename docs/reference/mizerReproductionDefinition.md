# Reuse the built-in reproduction capability definition

Returns the data-only balancing and reproduction-level capabilities for
the built-in law named in `params$rates_funcs$RDD`. This helper does not
register support or change the model. An extension author can explicitly
adopt the definition with
[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md)
after verifying that its effective reproduction route still implements
that law and that density-independent reproduction remains proportional
to `erepro`.

## Usage

``` r
mizerReproductionDefinition(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

A named list containing stationary, equation_rates_fun, tuning and
level_fun, suitable as the definition argument for the reproduction
process. The strategies reuse mizer's law-specific implementations.

## Details

The recognised laws are Beverton-Holt, Ricker, Shepherd, constant
recruitment and density-independent recruitment. Custom mortality or
encounter providers may leave these reproduction relationships
unchanged. A custom RDI, RDD or complete-rate wrapper needs its own
analysis: merely naming a built-in RDD does not establish that the
effective route preserves those relationships. Unknown law names are
rejected. Transformation support must be declared separately by the
extension author.

## See also

[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md),
[`balanceRDD()`](https://sizespectrum.org/mizer/reference/balanceRDD.md)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`other_mort()`](https://sizespectrum.org/mizer/reference/other_mort.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
