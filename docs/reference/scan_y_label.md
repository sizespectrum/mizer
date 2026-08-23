# Axis labels for a MizerScan

Assembles `"<name> [<units>]"`, the same way `array_y_label()` does for
the array classes.

## Usage

``` r
scan_y_label(x, default = "Value")

scan_x_label(x, default = "Scan value")

label_with_units(name, units)
```

## Arguments

- x:

  A MizerScan object.

- default:

  The label to use when the name is missing.

- name:

  The name of the quantity.

- units:

  The units, possibly NULL.

## Value

A string.
