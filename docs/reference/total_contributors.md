# Assemble the contributors to the total of a species-by-size array

The total is the total of everything the array holds: every species,
whether or not it was selected for display, and every size, whether or
not it falls in a species' own size range. It is a property of the array
rather than of the plot, so that a plot of two species can still be read
against the community total.

## Usage

``` r
total_contributors(x, wlim = c(NA, NA))
```

## Arguments

- x:

  An `ArraySpeciesBySize` object.

- wlim:

  Numeric vector of length two giving the weight limits.

## Value

A data frame of plotting data holding every value in the array.

## Details

The rows are returned unsummed, because the sum has to be taken after
the size coordinate has been converted — on a length axis the species no
longer share a grid; see
[`add_total_line()`](https://sizespectrum.org/mizer/reference/add_total_line.md).
