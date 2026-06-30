# Back-compatible wrapper for the core MizerParams upgrade

Retained because it was an (unexported) internal entry point. New code
should rely on
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
/
[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md),
which orchestrate the core upgrade together with any extension upgrades.

## Usage

``` r
upgradeParams(params)
```

## Arguments

- params:

  An old MizerParams object.

## Value

The upgraded MizerParams object.
