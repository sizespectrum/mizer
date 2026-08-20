# Factor relating a density measure to a density with respect to weight

Writing \\N_w\\ for the density with respect to weight, the density with
respect to measure \\m\\ is \\N_w\\ times the factor returned here. With
the allometric weight-length relationship \\w = a l^b\\ these factors
are \\1\\ for `"w"`, \\w\\ for `"log_w"`, \\dw/dl = b w / l\\ for `"l"`
and \\l\\dw/dl = b w\\ for `"log_l"`.

## Usage

``` r
density_measure_weight(measure, w, l, b)
```

## Arguments

- measure:

  One of
  [density_measures](https://sizespectrum.org/mizer/reference/density_measures.md).

- w, l, b:

  Numeric vectors of the same length giving the weight, the
  corresponding length, and the exponent of the weight-length
  relationship.

## Value

A numeric vector of factors.
