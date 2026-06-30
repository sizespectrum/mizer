# Sheperd function to calculate density-dependent reproduction rate

**\[experimental\]** Takes the density-independent rates \\R\_{di}\\ of
egg production and returns reduced, density-dependent rates \\R\_{dd}\\
given as \$\$R\_{dd} = \frac{R\_{di}}{1+(b\\ R\_{di})^c}\$\$

## Usage

``` r
SheperdRDD(rdi, species_params, ...)
```

## Arguments

- rdi:

  Vector of density-independent reproduction rates \\R\_{di}\\ for all
  species.

- species_params:

  A species parameter dataframe. Must contain columns `sheperd_b` and
  `sheperd_c` with the parameters b and c.

- ...:

  Unused

## Value

Vector of density-dependent reproduction rates.

## Details

With \\b = 1/R\_{max}\\ and \\c = 1\\ this reduces to the Beverton-Holt
reproduction rate, see
[`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md).

## See also

Other functions calculating density-dependent reproduction rate:
[`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md),
[`RickerRDD()`](https://sizespectrum.org/mizer/reference/RickerRDD.md),
[`constantEggRDI()`](https://sizespectrum.org/mizer/reference/constantEggRDI.md),
[`constantRDD()`](https://sizespectrum.org/mizer/reference/constantRDD.md),
[`noRDD()`](https://sizespectrum.org/mizer/reference/noRDD.md)
