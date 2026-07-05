# Length based sigmoid selectivity function

A sigmoid shaped selectivity function. Based on two parameters `l25` and
`l50` which determine the length at which 25% and 50% of the stock is
selected respectively.

## Usage

``` r
sigmoid_length(w, l25, l50, species_params, ...)
```

## Arguments

- w:

  Vector of sizes.

- l25:

  the length which gives a selectivity of 25%.

- l50:

  the length which gives a selectivity of 50%.

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
to `"sigmoid_length"` and provide the `l25` and `l50` values as
additional columns.
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
will then call this function automatically when calculating the
selectivity array.

The selectivity is given by the logistic function \$\$S(l) =
\frac{1}{1 + \exp\left(\log(3)\frac{l50 -l}{l50 - l25}\right)}\$\$ As
the mizer model is weight based, and this selectivity function is length
based, it uses the length-weight parameters `a` and `b` to convert
between length and weight \$\$l = \left(\frac{w}{a}\right)^{1/b}\$\$

## See also

[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
for setting the selectivity parameters.

Other selectivity functions:
[`double_sigmoid_length()`](https://sizespectrum.org/mizer/reference/double_sigmoid_length.md),
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md),
[`sigmoid_weight()`](https://sizespectrum.org/mizer/reference/sigmoid_weight.md)

## Examples

``` r
# Selectivity at weight given l25 = 10 cm, l50 = 15 cm
# using length-weight parameters a = 0.01, b = 3
sp <- list(a = 0.01, b = 3)
w <- c(1, 10, 100, 500, 1000)
sigmoid_length(w, l25 = 10, l50 = 15, species_params = sp)
#> [1] 0.0931323 0.2500000 0.8081354 0.9918278 0.9989960
```
