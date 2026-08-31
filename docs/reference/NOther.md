# Time series of other components

Fetch the simulation results for other components over time.

## Usage

``` r
NOther(sim)

finalNOther(sim)
```

## Arguments

- sim:

  A MizerSim object

## Value

For `NOther`: A list array indexed by time and component that stores the
projected values for other ecosystem components.

For `finalNOther`: A named list holding the values of other ecosystem
components at the end of the simulation

## See also

Other extension tools:
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`other_mort()`](https://sizespectrum.org/mizer/reference/other_mort.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
