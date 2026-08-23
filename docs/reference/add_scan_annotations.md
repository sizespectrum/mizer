# Add the annotation layers to a MizerScan plot

Add the annotation layers to a MizerScan plot

## Usage

``` r
add_scan_annotations(
  p,
  x,
  plot_dat,
  reference_lines = TRUE,
  mark_max = FALSE,
  show_unsettled = TRUE
)
```

## Arguments

- p:

  The ggplot object so far.

- x:

  The MizerScan object.

- plot_dat:

  The data frame being plotted.

- reference_lines:

  TRUE to use the stored reference lines, FALSE for none, or a named
  numeric vector to use instead.

- mark_max:

  Whether to mark where each series attains its maximum.

- show_unsettled:

  Whether to mark the scan values that did not settle.

## Value

The ggplot object with the extra layers, still a `mizer_plot`.
