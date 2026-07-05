# Beverton Holt function to calculate density-dependent reproduction rate

Takes the density-independent rates \\R\_{di}\\ of egg production (as
calculated by
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md)) and
returns reduced, density-dependent reproduction rates \\R\_{dd}\\ given
as \$\$R\_{dd} = R\_{di} \frac{R\_{max}}{R\_{di} + R\_{max}}\$\$ where
\\R\_{max}\\ are the maximum possible reproduction rates that must be
specified in a column in the species parameter dataframe. (All
quantities in the above equation are species-specific but we dropped the
species index for simplicity.)

## Usage

``` r
BevertonHoltRDD(rdi, species_params, ...)
```

## Arguments

- rdi:

  Vector of density-independent reproduction rates \\R\_{di}\\ for all
  species.

- species_params:

  A species parameter dataframe. Must contain a column `R_max` holding
  the maximum reproduction rate \\R\_{max}\\ for each species.

- ...:

  Unused

## Value

Vector of density-dependent reproduction rates.

## Details

This is only one example of a density-dependence. You can write your own
function based on this example, returning different density-dependent
reproduction rates. Three other examples provided are
[`RickerRDD()`](https://sizespectrum.org/mizer/reference/RickerRDD.md),
[`SheperdRDD()`](https://sizespectrum.org/mizer/reference/SheperdRDD.md),
[`noRDD()`](https://sizespectrum.org/mizer/reference/noRDD.md) and
[`constantRDD()`](https://sizespectrum.org/mizer/reference/constantRDD.md).
For more explanation see
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md).

## See also

Other functions calculating density-dependent reproduction rate:
[`RickerRDD()`](https://sizespectrum.org/mizer/reference/RickerRDD.md),
[`SheperdRDD()`](https://sizespectrum.org/mizer/reference/SheperdRDD.md),
[`constantEggRDI()`](https://sizespectrum.org/mizer/reference/constantEggRDI.md),
[`constantRDD()`](https://sizespectrum.org/mizer/reference/constantRDD.md),
[`noRDD()`](https://sizespectrum.org/mizer/reference/noRDD.md)
