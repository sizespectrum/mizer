# Measure distance between current and previous state in terms of fish abundances

**\[experimental\]**

Calculates the sum squared difference between log(N) in current and
previous state. This function can be used in
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
to decide when sufficient convergence to steady state has been achieved.

Only the size classes that hold fish are measured. A class whose density
is zero in either state has no log to take, and one holding a negligible
share of its species' biomass has a log that never stops moving: above a
size where growth stops the density decays exponentially, so `log(n)`
falls by the same amount between every pair of states, for ever. Left in
the sum, a single such class holding \\10^{-92}\\ g of fish can hold the
distance above any tolerance indefinitely and stop
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
from ever converging, while the biomass drift correctly reports a fixed
point. Excluding it changes nothing for a model that does not have one,
because the classes that hold the fish are exactly the classes that are
kept.

## Usage

``` r
distanceSSLogN(
  params,
  current,
  previous,
  biomass_share_cutoff = steady_share_cutoff(),
  ...
)
```

## Arguments

- params:

  MizerParams

- current:

  A named list with entries `n`, `n_pp` and `n_other` describing the
  current state

- previous:

  A named list with entries `n`, `n_pp` and `n_other` describing the
  previous state

- biomass_share_cutoff:

  **\[experimental\]** A finite number between 0 and 1, inclusive,
  giving the share of a species' biomass that a size class must hold in
  the current state to be measured. `0` measures every class with a
  nonzero density in both states, which is what this function did before
  mizer 3.3.

- ...:

  Unused. Accepted because
  [`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
  forwards its own `...` to both the rate functions and the distance
  function, so an argument meant for one arrives at the other.

## Value

The sum of squares of the difference in the logs of the fish abundances
`n`, over the size classes that hold fish in both states:
`sum((log(current$n) - log(previous$n))^2)`

## See also

Other distance functions:
[`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md)
