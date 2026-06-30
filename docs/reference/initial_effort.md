# Initial fishing effort

The fishing effort is a named vector, specifying for each fishing gear
the effort invested into fishing with that gear. The effort value for
each gear is multiplied by the catchability and the selectivity to
determine the fishing mortality imposed by that gear, see
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
for more details. The initial effort you have set can be overruled when
running a simulation by providing an `effort` argument to
[`project()`](https://sizespectrum.org/mizer/reference/project.md) which
allows you to specify a time-varying effort.

## Usage

``` r
initial_effort(params)

initial_effort(params) <- value
```

## Arguments

- params:

  A MizerParams object

- value:

  A vector or scalar with the initial fishing effort, see Details below.

## Value

A named effort vector ordered by gear.

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
  by
  [`validEffortVector()`](https://sizespectrum.org/mizer/reference/validEffortVector.md),
  which depends on the defaults edition.

These conversions are done by the function
[`validEffortVector()`](https://sizespectrum.org/mizer/reference/validEffortVector.md).

An `effort` argument will lead to an error if it is either

- unnamed and of the wrong length

- named but where some names do not match any of the gears

- not numeric
