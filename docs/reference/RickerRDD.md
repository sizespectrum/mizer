# Ricker function to calculate density-dependent reproduction rate

**\[experimental\]** Takes the density-independent rates \\R\_{di}\\ of
egg production and returns reduced, density-dependent rates \\R\_{dd}\\
given as \$\$R\_{dd} = R\_{di} \exp(- b R\_{di})\$\$

## Usage

``` r
RickerRDD(rdi, species_params, ...)
```

## Arguments

- rdi:

  Vector of density-independent reproduction rates \\R\_{di}\\ for all
  species.

- species_params:

  A species parameter dataframe. Must contain a column `ricker_b`
  holding the coefficient b.

- ...:

  Unused

## Value

Vector of density-dependent reproduction rates.

## See also

Other functions calculating density-dependent reproduction rate:
[`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md),
[`SheperdRDD()`](https://sizespectrum.org/mizer/reference/SheperdRDD.md),
[`constantEggRDI()`](https://sizespectrum.org/mizer/reference/constantEggRDI.md),
[`constantRDD()`](https://sizespectrum.org/mizer/reference/constantRDD.md),
[`noRDD()`](https://sizespectrum.org/mizer/reference/noRDD.md)
