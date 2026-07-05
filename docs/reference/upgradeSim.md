# Back-compatible wrapper for the MizerSim upgrade

Retained because it was an (unexported) internal entry point. New code
should rely on
[`validSim()`](https://sizespectrum.org/mizer/reference/validSim.md) /
[`readSim()`](https://sizespectrum.org/mizer/reference/saveParams.md).

## Usage

``` r
upgradeSim(sim)
```

## Arguments

- sim:

  An old MizerSim object.

## Value

The upgraded MizerSim object.
