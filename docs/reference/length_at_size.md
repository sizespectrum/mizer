# Length of an individual at each weight on the size grid

Internal helper for
[`getMeanLength()`](https://sizespectrum.org/mizer/reference/getMeanWeight.md).
Inverts the length-weight relationship \\w = a_i l^{b_i}\\ of each
species at every weight in `params@w`. The result is a weighting factor
for
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md),
so it is a point value on the grid and is not bin-averaged here.

## Usage

``` r
length_at_size(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

A species x size matrix of lengths in cm.
