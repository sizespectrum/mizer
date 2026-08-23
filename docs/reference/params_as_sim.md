# A MizerSim holding only the state stored in a MizerParams

Wraps the single state in a `MizerParams` object into a `MizerSim` with
one time step, so that a function written against a `MizerSim` can be
applied to it without projecting.

## Usage

``` r
params_as_sim(params, t = 0)
```

## Arguments

- params:

  A MizerParams object.

- t:

  The time to label the single time step with.

## Value

A MizerSim object with one time step holding the initial state of
`params`.

## Details

This exists for analyses such as
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
that measure a quantity with a user-supplied function of a `MizerSim`.
When the model has settled on a fixed point there is nothing to project:
the state does not change, so a snapshot of it carries all the
information a longer run would. Every slot that a summary function might
read has to be filled, because
[`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md)
initialises them all to `NA` and, for example,
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md)
reads `sim@effort` and would otherwise return `NA` without complaint.

Note that the result has a single time step, so a function that needs
more than one — anything taking a difference through time — cannot be
applied to it.
