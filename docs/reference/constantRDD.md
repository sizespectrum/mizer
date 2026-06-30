# Give constant reproduction rate

**\[experimental\]** Simply returns the value from
`species_params$constant_reproduction`.

## Usage

``` r
constantRDD(rdi, species_params, ...)
```

## Arguments

- rdi:

  Vector of density-independent reproduction rates \\R\_{di}\\ for all
  species.

- species_params:

  A species parameter dataframe. Must contain a column
  `constant_reproduction`.

- ...:

  Unused

## Value

Vector `species_params$constant_reproduction`

## See also

Other functions calculating density-dependent reproduction rate:
[`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md),
[`RickerRDD()`](https://sizespectrum.org/mizer/reference/RickerRDD.md),
[`SheperdRDD()`](https://sizespectrum.org/mizer/reference/SheperdRDD.md),
[`constantEggRDI()`](https://sizespectrum.org/mizer/reference/constantEggRDI.md),
[`noRDD()`](https://sizespectrum.org/mizer/reference/noRDD.md)
