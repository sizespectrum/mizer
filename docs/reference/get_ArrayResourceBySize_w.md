# Get the size grid for an ArrayResourceBySize object

Internal helper that follows the stored mapping to the full
prey/resource size grid. Falls back to numeric names when no model grid
is available.

## Usage

``` r
get_ArrayResourceBySize_w(x)
```

## Arguments

- x:

  An `ArrayResourceBySize` object.

## Value

A numeric vector giving the size represented by each element.
