# Strip extension classes from a mizer object

Coerces a `MizerParams` or `MizerSim` object back to its plain base
class, removing any S3 extension classes. For `MizerSim`, also strips
the extension class from the embedded `params` element.

## Usage

``` r
baseMizerClass(object)
```

## Arguments

- object:

  A `MizerParams` or `MizerSim` object.

## Value

The same object coerced to `MizerParams` or `MizerSim`.
