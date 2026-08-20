# Length based knife-edge selectivity function

A knife-edge selectivity function where individuals with a length
greater or equal to `knife_edge_length` are fully selected and no fish
shorter than this length are selected.

## Usage

``` r
knife_edge_length(w, knife_edge_length, species_params, ...)
```

## Arguments

- w:

  Vector of sizes (weights).

- knife_edge_length:

  The length at which the knife-edge operates.

- species_params:

  A list with the species params for the current species. Used to get at
  the length-weight parameters `a` and `b`.

- ...:

  Unused

## Value

Vector of selectivities at the given sizes.

## Details

You would not usually call this function directly. Instead, set the
`sel_func` column in
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
to `"knife_edge_length"` and provide `knife_edge_length` as an
additional column.
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
will then call this function automatically when calculating the
selectivity array.

As the mizer model is weight based, and this selectivity function is
length based, it uses the length-weight parameters `a` and `b` to
convert the cut-off length to a weight: \$\$w\_{\text{cut}} = a \cdot
l\_{\text{cut}}^b\$\$

## See also

[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
for setting the `knife_edge_length` parameter.

Other selectivity functions:
[`double_sigmoid_length()`](https://sizespectrum.org/mizer/reference/double_sigmoid_length.md),
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md),
[`sigmoid_length()`](https://sizespectrum.org/mizer/reference/sigmoid_length.md),
[`sigmoid_weight()`](https://sizespectrum.org/mizer/reference/sigmoid_weight.md)

## Examples

``` r
# Knife-edge at 20 cm using length-weight parameters a = 0.01, b = 3
sp <- list(a = 0.01, b = 3)
knife_edge_length(w = c(1, 10, 100, 1000), knife_edge_length = 20,
                  species_params = sp)
#> [1] 0 0 1 1
```
