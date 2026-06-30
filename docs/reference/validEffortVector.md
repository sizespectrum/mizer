# Make a valid effort vector

Make a valid effort vector

## Usage

``` r
validEffortVector(effort, params)
```

## Arguments

- effort:

  A vector or scalar with the initial fishing effort, see Details below.

- params:

  A MizerParams object.

## Details

A valid effort vector is a named vector with one effort value for each
gear. However you can also supply the effort value in different ways:

- a scalar, which is then replicated for each gear

- an unnamed vector, which is then assumed to be in the same order as
  the gears in the params object

- a named vector in which the gear names have a different order than in
  the params object. This is then sorted correctly.

- a named vector which only supplies values for some of the gears. The
  effort for the other gears is then set to the default effort returned
  by `validEffortVector()`, which depends on the defaults edition.

These conversions are done by the function `validEffortVector()`.

An `effort` argument will lead to an error if it is either

- unnamed and of the wrong length

- named but where some names do not match any of the gears

- not numeric

## See also

[`initial_effort()`](https://sizespectrum.org/mizer/reference/initial_effort.md)
