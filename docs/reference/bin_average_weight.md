# Trapezoidal bin-average of a per-bin weight

Internal helper for the second-order summary integrals. A summary
diagnostic \\\int N(w) K(w)\\ dw\\ is discretised on the finite-volume
grid as \\\sum_j N_j \bar K_j \Delta w_j\\, where \\N_j\\ is the cell
average of the density over bin \\\[w_j, w\_{j+1}\]\\. To be second
order in the bin width the point weight \\K(w_j)\\ must be replaced by
the bin average \$\$\bar K_j = \frac{1}{\Delta
w_j}\int\_{w_j}^{w\_{j+1}} K(w)\\dw \approx \tfrac12\big(K(w_j) +
K(w\_{j+1})\big).\$\$ The trapezoidal average \\\tfrac12(K_j +
K\_{j+1})\\ is uniformly second order and exact whenever \\K\\ is linear
in \\w\\ (e.g. the first moment \\K = w\\, for which it equals
\\(w\_{j+1}^2 - w_j^2)/(2\Delta w_j)\\).

## Usage

``` r
bin_average_weight(K)
```

## Arguments

- K:

  A numeric vector of weights indexed over the size grid, or a numeric
  array whose last dimension runs over the size grid (e.g. a
  species-by-size matrix or a gear-by-species-by-size array).

## Value

An object of the same shape as `K` containing the trapezoidal
bin-averaged weights.

## Details

The weight `K` is supplied already evaluated on the size grid (a vector
indexed over the bins, or a matrix with the size dimension running along
the columns). The top bin has no right-hand neighbour on the grid, so
its weight is left unaveraged (one-sided); the density there is
negligible, so this does not affect the second-order accuracy of the
totals.

This helper is shared with the reproduction integrals (issue \#376),
which also need the trapezoidal bin-average of a composite weight.
