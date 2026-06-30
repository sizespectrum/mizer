# Measure distance between current and previous state in terms of fish abundances

**\[experimental\]**

Calculates the sum squared difference between log(N) in current and
previous state. This function can be used in
[`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md)
to decide when sufficient convergence to steady state has been achieved.

## Usage

``` r
distanceSSLogN(params, current, previous)
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

## Value

The sum of squares of the difference in the logs of the (nonzero) fish
abundances `n`, ignoring entries where either state has zero abundance:
`sum((log(current$n) - log(previous$n))^2)`

## See also

Other distance functions:
[`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md)
