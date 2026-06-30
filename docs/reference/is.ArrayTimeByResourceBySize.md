# Test if an object is an ArrayTimeByResourceBySize

Test if an object is an ArrayTimeByResourceBySize

## Usage

``` r
is.ArrayTimeByResourceBySize(x)
```

## Arguments

- x:

  An object to test.

## Value

`TRUE` if `x` is an `ArrayTimeByResourceBySize` object, `FALSE`
otherwise.

## Examples

``` r
is.ArrayTimeByResourceBySize(NResource(NS_sim))
#> [1] TRUE
is.ArrayTimeByResourceBySize(matrix(1:4, nrow = 2))
#> [1] FALSE
```
