# Determine the tooltip variables for a mizer plot

Works out which variables should appear in the plotly tooltip, including
the legend variable only when it differs from the grouping variable,
plus any extra variables.

## Usage

``` r
mizer_tooltip_vars(
  frame,
  group_var,
  x_var,
  y_var,
  legend_var = NULL,
  extra = NULL
)
```

## Arguments

- frame:

  The data frame underlying the plot.

- group_var:

  Name of the grouping variable.

- x_var:

  Name of the variable on the x-axis.

- y_var:

  Name of the variable on the y-axis.

- legend_var:

  Optional name of the legend variable.

- extra:

  Optional character vector of additional variable names to include.

## Value

A character vector of variable names for the tooltip.
