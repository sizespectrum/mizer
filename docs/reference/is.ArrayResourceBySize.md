# Test if an object is an ArrayResourceBySize

Test if an object is an ArrayResourceBySize

## Usage

``` r
is.ArrayResourceBySize(x)
```

## Arguments

- x:

  An object to test.

## Value

`TRUE` if `x` is an `ArrayResourceBySize` object, `FALSE` otherwise.

## Examples

``` r
is.ArrayResourceBySize(getResourceMort(NS_params))
#> [1] TRUE
is.ArrayResourceBySize(1:4)
#> [1] FALSE
```
