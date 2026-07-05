# Test if an object is a ArrayTimeBySpecies

Test if an object is a ArrayTimeBySpecies

## Usage

``` r
is.ArrayTimeBySpecies(x)
```

## Arguments

- x:

  An object to test.

## Value

`TRUE` if `x` is an `ArrayTimeBySpecies` object, `FALSE` otherwise.

## Examples

``` r
is.ArrayTimeBySpecies(getBiomass(NS_sim))
#> [1] TRUE
is.ArrayTimeBySpecies(matrix(1:4, nrow = 2))
#> [1] FALSE
```
