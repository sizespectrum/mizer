# Prepare the data frame for plotting a MizerScan

Prepare the data frame for plotting a MizerScan

## Usage

``` r
prepare_MizerScan_plot_data(x, species = NULL)
```

## Arguments

- x:

  A MizerScan object.

- species:

  The series to keep, or NULL for all of them.

## Value

A data frame with the x, y and grouping variable in the first three
columns, as
[`plotDataFrame()`](https://sizespectrum.org/mizer/reference/plotDataFrame.md)
requires.
