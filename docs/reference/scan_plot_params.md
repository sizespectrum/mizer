# The params object to use when plotting a MizerScan

Series that are not species have no colour in the model, and
[`plotDataFrame()`](https://sizespectrum.org/mizer/reference/plotDataFrame.md)
silently drops any legend level it cannot find a colour for. So any such
series is given a colour here, using the ordinary
[`setColours()`](https://sizespectrum.org/mizer/reference/setColours.md)
interface, which also leaves the user free to choose a different one.

## Usage

``` r
scan_plot_params(x, plot_dat)
```

## Arguments

- x:

  A MizerScan object.

- plot_dat:

  The data frame that will be plotted.

## Value

A MizerParams object with a colour for every series in `plot_dat`.
