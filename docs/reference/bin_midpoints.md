# Geometric bin centres of the size grid

Internal helper for the second-order plotting code. A finite-volume cell
average \\N_j = (1/\Delta w_j)\int\_{w_j}^{w\_{j+1}} N\\dw\\ does not
live at the left bin edge \\w_j\\ but at the geometric bin centre
\$\$w^\*\_j = \sqrt{w_j\\w\_{j+1}} = w_j\sqrt{\beta},\$\$ where \\\beta
= w\_{j+1}/w_j\\ is the (constant) bin ratio of the logarithmic grid.
This is the log-symmetric, second-order-correct location at which to
plot a bin-averaged quantity (it is exact for the community spectrum
\\N\propto w^{-2}\\). It is a uniform half-bin shift to the right on the
log axis, the same for the consumer grid `w` and the full prey/resource
grid `w_full`.

## Usage

``` r
bin_midpoints(params, w_full = FALSE)
```

## Arguments

- params:

  A MizerParams object.

- w_full:

  If `TRUE`, return the centres of the full (resource) grid
  `params@w_full`; otherwise the consumer grid `params@w`.

## Value

A numeric vector of geometric bin centres, one per grid node.
