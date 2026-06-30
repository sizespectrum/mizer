# Test if an object is a ArraySpeciesBySize

Test if an object is a ArraySpeciesBySize

## Usage

``` r
is.ArraySpeciesBySize(x)
```

## Arguments

- x:

  An object to test.

## Value

`TRUE` if `x` is an `ArraySpeciesBySize` object, `FALSE` otherwise.

## Examples

``` r
is.ArraySpeciesBySize(getEncounter(NS_params))
#> [1] TRUE
is.ArraySpeciesBySize(matrix(1:4, nrow = 2))
#> [1] FALSE
```
