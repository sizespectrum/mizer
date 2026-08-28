# Bin-average the weight of a size-spectrum integral

**\[experimental\]** Prepares the weight \\K(w)\\ of an integral over
the size spectrum so that the integral is evaluated with the quadrature
scheme the model is actually using. Use this when writing your own
indicator or diagnostic function; the built-in summary and indicator
functions call it for you.

## Usage

``` r
bin_average_weight(K, params)
```

## Arguments

- K:

  A numeric vector of weights indexed over the size grid, or a numeric
  array whose last dimension runs over the size grid (e.g. a
  species-by-size matrix or a gear-by-species-by-size array).

- params:

  A MizerParams object whose `second_order_w` slot controls the gating.

## Value

The weight `K`, bin-averaged when
`params@second_order_w[["bin_average"]]` is `TRUE`, otherwise returned
unchanged.

## Details

An integral \\\int N(w) K(w)\\ dw\\ is discretised on mizer's
finite-volume grid as \\\sum_j N_j \bar K_j \Delta w_j\\, where \\N_j\\
is the cell average of the density over bin \\\[w_j, w\_{j+1}\]\\. Only
the weight is approximated: \\N_j\\ is already a cell average and
\\\Delta w_j\\ is exact, so **neither the abundance nor the bin widths
should ever be passed through this function**.

Whether the point weight \\K(w_j)\\ is replaced by the bin average
\$\$\bar K_j = \frac{1}{\Delta w_j}\int\_{w_j}^{w\_{j+1}} K(w)\\dw
\approx \tfrac12\big(K(w_j) + K(w\_{j+1})\big)\$\$ is controlled by the
`bin_average` entry of the model's
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
slot. When it is `FALSE` (the default) `K` is returned unchanged, so an
indicator written with this function reproduces the left-edge Riemann
sums of previous mizer versions byte-for-byte. When it is `TRUE` the
trapezoidal bin average is returned, which is uniformly second order and
exact whenever \\K\\ is linear in \\w\\ (e.g. the first moment \\K =
w\\, for which it equals \\(w\_{j+1}^2 - w_j^2)/(2\Delta w_j)\\).

Because the gating happens inside, always call this rather than
averaging unconditionally: a hard-coded bin average silently changes the
results of models that are on the default scheme.

If `K` is a product of several size-dependent factors, average the
**product** and not the individual factors — the average of a product is
not the product of the averages. Spawning stock biomass, for example,
averages `maturity * w` as a single weight.

The top bin has no right-hand neighbour on the grid, so its weight is
left unaveraged (one-sided); the density there is negligible, so this
does not affect the second-order accuracy of the totals.

## See also

[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md),
[`get_size_range_array()`](https://sizespectrum.org/mizer/reference/get_size_range_array.md),
[`encounter_kernel()`](https://sizespectrum.org/mizer/reference/encounter_kernel.md)

## Examples

``` r
# Biomass of each species above 10g -- what getBiomass() does internally.
params <- NS_params
K <- get_size_range_array(params, min_w = 10)   # species x size, 0/1
K <- sweep(K, 2, params@w, "*")                 # weight by w to get biomass
K <- bin_average_weight(K, params)              # gated on second_order_w
rowSums(sweep(initialN(params) * K, 2, params@dw, "*"))
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#>        Sprat      Sandeel       N.pout      Herring          Dab      Whiting 
#> 2.441192e+11 4.589606e+12 2.388476e+11 1.273446e+12 8.373096e+09 1.622396e+11 
#>         Sole      Gurnard       Plaice      Haddock          Cod       Saithe 
#> 1.274005e+11 2.475989e+10 7.694543e+11 3.559899e+11 5.998285e+11 4.890086e+11 
```
