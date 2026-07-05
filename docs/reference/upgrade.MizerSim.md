# Upgrade a MizerSim object from earlier versions

This is the `MizerSim` method of the
[`upgrade()`](https://rdrr.io/r/utils/upgrade.html) generic. It rebuilds
the simulation around an upgraded params object; the params upgrade
(core mizer and any extensions) is performed by the
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
call below. You should never need to call it directly; use
[`validSim()`](https://sizespectrum.org/mizer/reference/validSim.md) (or
[`readSim()`](https://sizespectrum.org/mizer/reference/saveParams.md)).

## Usage

``` r
# S3 method for class 'MizerSim'
upgrade(object, ...)
```

## Arguments

- object:

  An old MizerSim object to be upgraded

- ...:

  Unused.

## Value

The upgraded MizerSim object
