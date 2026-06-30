# Weight based knife-edge selectivity function

A knife-edge selectivity function where weights greater or equal to
`knife_edge_size` are fully selected and no fish smaller than this size
are selected.

## Usage

``` r
knife_edge(w, knife_edge_size, ...)
```

## Arguments

- w:

  Vector of sizes.

- knife_edge_size:

  The weight at which the knife-edge operates.

- ...:

  Unused

## Value

Vector of selectivities at the given sizes.

## Details

You would not usually call this function directly. Instead, set the
`sel_func` column in
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
to `"knife_edge"` and provide `knife_edge_size` as an additional column.
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
will then call this function automatically when calculating the
selectivity array.

## See also

[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
for setting the `knife_edge_size` parameter.

Other selectivity functions:
[`double_sigmoid_length()`](https://sizespectrum.org/mizer/reference/double_sigmoid_length.md),
[`sigmoid_length()`](https://sizespectrum.org/mizer/reference/sigmoid_length.md),
[`sigmoid_weight()`](https://sizespectrum.org/mizer/reference/sigmoid_weight.md)

## Examples

``` r
knife_edge(w = c(1, 10, 100, 1000), knife_edge_size = 100)
#> [1] 0 0 1 1
```
