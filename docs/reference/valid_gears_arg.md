# Helper function to assure validity of gears argument

If the gears argument contains invalid gears, then these are ignored but
a warning is issued.

## Usage

``` r
valid_gears_arg(object, gears = NULL, error_on_empty = FALSE)
```

## Arguments

- object:

  A MizerSim or MizerParams object from which the gears should be
  selected.

- gears:

  The gears to be selected. Optional. By default all gears are selected.
  A vector of gear names.

- error_on_empty:

  Whether to throw an error if there are zero valid gears. Default
  FALSE.

## Value

A vector of gear names in the same order as supplied in `gears`, with
invalid names removed. If `gears` is `NULL`, all gears are returned in
the order stored in the model.
