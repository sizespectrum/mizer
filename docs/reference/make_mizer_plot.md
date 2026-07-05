# Tag a ggplot object as a mizer plot

Attaches the tooltip information to a ggplot object and adds the
`"mizer_plot"` class so that
[`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md)
knows which aesthetics to show.

## Usage

``` r
make_mizer_plot(plot, tooltip)
```

## Arguments

- plot:

  A ggplot object.

- tooltip:

  Character vector of variable names to include in the plotly tooltip.

## Value

The `plot` object with the tooltip stored as an attribute and
`"mizer_plot"` prepended to its class.
