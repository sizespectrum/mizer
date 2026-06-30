# Weight based sigmoidal selectivity function

A sigmoidal selectivity function with 50% selectivity at weight
`sigmoidal_weight` \\=w\_{\text{sigmoid}}\\ and width `sigmoidal_sigma`
\\=\sigma\\. \$\$S(w) = \left(1 +
\left(\frac{w}{w\_{\text{sigmoid}}}\right)^{-\sigma}\right)^{-1}\$\$

## Usage

``` r
sigmoid_weight(w, sigmoidal_weight, sigmoidal_sigma, ...)
```

## Arguments

- w:

  Vector of sizes.

- sigmoidal_weight:

  The weight at which selectivity is 50%.

- sigmoidal_sigma:

  The width of the selection function.

- ...:

  Unused

## Value

Vector of selectivities at the given sizes.

## Details

You would not usually call this function directly. Instead, set the
`sel_func` column in
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
to `"sigmoid_weight"` and provide `sigmoidal_weight` and
`sigmoidal_sigma` as additional columns.
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
will then call this function automatically when calculating the
selectivity array.

## See also

[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
for setting the selectivity parameters.

Other selectivity functions:
[`double_sigmoid_length()`](https://sizespectrum.org/mizer/reference/double_sigmoid_length.md),
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md),
[`sigmoid_length()`](https://sizespectrum.org/mizer/reference/sigmoid_length.md)

## Examples

``` r
sigmoid_weight(w = c(1, 10, 100, 1000),
               sigmoidal_weight = 100, sigmoidal_sigma = 3)
#> [1] 9.99999e-07 9.99001e-04 5.00000e-01 9.99001e-01
```
