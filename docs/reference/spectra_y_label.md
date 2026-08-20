# Y-axis label for a size-spectrum plot

Y-axis label for a size-spectrum plot

## Usage

``` r
spectra_y_label(
  power,
  size_axis = "w",
  biomass = power >= 1,
  per_log_size = power == 2
)
```

## Arguments

- power:

  The power of weight that the abundance was multiplied by.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- biomass:

  Whether the quantity is a biomass density rather than a number
  density. Defaults to `power >= 1`.

- per_log_size:

  Whether the quantity is a density with respect to logarithmic size.
  Defaults to `power == 2`.

## Value

A character string for the y-axis label.
