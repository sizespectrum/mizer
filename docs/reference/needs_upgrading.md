# Determine whether a MizerParams or MizerSim object needs to be upgraded

Looks at the mizer version that was used to last update the object and
returns TRUE if changes since that version require an upgrade of the
object. You would not usually have to call this function. Upgrades are
initiated automatically by `validParams` and `validSim` when necessary.

## Usage

``` r
needs_upgrading(object)
```

## Arguments

- object:

  A MizerParams or MizerSim object

## Value

TRUE or FALSE
