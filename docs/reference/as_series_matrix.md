# Bring what a value function returned into a time by series matrix

Bring what a value function returned into a time by series matrix

## Usage

``` r
as_series_matrix(x, default_name = "Value")
```

## Arguments

- x:

  The return value of a `value_func`.

- default_name:

  The column name to use for a single unnamed series.

## Value

A matrix with time in the rows and the series in the columns.
