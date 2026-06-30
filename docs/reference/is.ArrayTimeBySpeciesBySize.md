# Test if an object is an ArrayTimeBySpeciesBySize

Test if an object is an ArrayTimeBySpeciesBySize

## Usage

``` r
is.ArrayTimeBySpeciesBySize(x)
```

## Arguments

- x:

  An object to test.

## Value

`TRUE` if `x` is an `ArrayTimeBySpeciesBySize` object, `FALSE`
otherwise.

## Examples

``` r
is.ArrayTimeBySpeciesBySize(getFMort(NS_sim))
#> [1] TRUE
is.ArrayTimeBySpeciesBySize(array(1:8, dim = c(2, 2, 2)))
#> [1] FALSE
```
