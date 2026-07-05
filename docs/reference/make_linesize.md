# Construct a named vector of line widths for a plot

Helper used by the plotting functions to give highlighted species a
thicker line than the rest.

## Usage

``` r
make_linesize(levels, highlight)
```

## Arguments

- levels:

  Character vector of the legend levels (usually species names).

- highlight:

  Name or vector of names of the levels to be highlighted with a thicker
  line.

## Value

A named numeric vector of line widths, one for each entry in `levels`,
with highlighted entries set to a larger value.
