# Whether a plot's y axis is logarithmic

Read from the plot's own y scale, so that data being *added* to an
existing plot can be filtered the way that plot was. A plot with no
explicit y scale has ggplot2's default, which is linear.

## Usage

``` r
plot_y_is_log(plot)
```

## Arguments

- plot:

  A ggplot2 object.

## Value

`TRUE` if the plot's y axis uses a log10 transformation.
