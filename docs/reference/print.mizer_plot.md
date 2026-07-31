# Print a mizer plot

Suppresses the uninformative ggplot2 warning about log transformations
introducing infinite values, which occurs when zero values are present
on a logged axis.

## Usage

``` r
# S3 method for class 'mizer_plot'
print(x, ...)
```

## Arguments

- x:

  A `mizer_plot` object.

- ...:

  Further arguments passed to the ggplot2 print method.

## Value

The plot object, invisibly.
