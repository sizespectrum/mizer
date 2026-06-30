# Run the registered extension upgrade methods on an object

For each extension recorded in the object's `@extensions` slot
(processed innermost-first) whose installed package version is newer
than the recorded stamp (or whose stamp is missing), calls the
extension's `upgrade` method if one is registered, then records the
installed version as the new stamp. The core mizer upgrade is *not* run
here; see
[`upgrade.MizerParams()`](https://sizespectrum.org/mizer/reference/upgrade.MizerParams.md).

## Usage

``` r
runExtensionUpgrades(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

The object with extension migrations applied and stamps refreshed.

## Details

Extension upgrade methods are looked up with
`getS3method("upgrade", name, optional = TRUE)`, must perform only their
own migration, must be idempotent, and must not call
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).
