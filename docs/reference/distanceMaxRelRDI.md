# Measure distance between current and previous state in terms of RDI

**\[experimental\]**

This function can be used in
[`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md)
to decide when sufficient convergence to steady state has been achieved.

## Usage

``` r
distanceMaxRelRDI(params, current, previous)
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

The largest absolute relative change in rdi:
`max(abs((current_rdi - previous_rdi) / previous_rdi))`. If any entry of
`previous_rdi` is zero, the result can be infinite.

## See also

Other distance functions:
[`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md)
