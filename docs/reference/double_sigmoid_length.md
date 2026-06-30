# Length based double-sigmoid selectivity function

A hump-shaped selectivity function with a sigmoidal rise and an
independent sigmoidal drop-off. This drop-off is what distinguishes this
from the function
[`sigmoid_length()`](https://sizespectrum.org/mizer/reference/sigmoid_length.md)
and it is intended to model the escape of large individuals from the
fishing gear.

## Usage

``` r
double_sigmoid_length(w, l25, l50, l50_right, l25_right, species_params, ...)
```

## Arguments

- w:

  Vector of sizes.

- l25:

  the length which gives a selectivity of 25%.

- l50:

  the length which gives a selectivity of 50%.

- l50_right:

  the length which gives a selectivity of 50%.

- l25_right:

  the length which gives a selectivity of 25%.

- species_params:

  A list with the species params for the current species. Used to get at
  the length-weight parameters `a` and `b`

- ...:

  Unused

## Value

Vector of selectivities at the given sizes. Requires
`l25 < l50 < l50_right < l25_right`.

## Details

You would not usually call this function directly. Instead, set the
`sel_func` column in
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
to `"double_sigmoid_length"` and provide the `l25`, `l50`, `l50_right`
and `l25_right` values as additional columns.
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
will then call this function automatically when calculating the
selectivity array.

The selectivity is obtained as the product of two sigmoidal curves, one
rising and one dropping. The sigmoidal rise is based on the two
parameters `l25` and `l50` which determine the length at which 25% and
50% of the stock is selected respectively. The sigmoidal drop-off is
based on the two parameters `l50_right` and `l25_right` which determine
the length at which the selectivity curve has dropped back to 50% and
25% respectively. The selectivity is given by the function \$\$S(l) =
\frac{1}{1 + \exp\left(\log(3)\frac{l50 -l}{l50 -
l25}\right)}\frac{1}{1 + \exp\left(\log(3)\frac{l50\_{right}
-l}{l50\_{right} - l25\_{right}}\right)}\$\$

As the size-based model is weight based, and this selectivity function
is length based, it uses the length-weight parameters `a` and `b` to
convert between length and weight. \$\$l =
\left(\frac{w}{a}\right)^{1/b}\$\$

## See also

[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
for setting the selectivity parameters.

Other selectivity functions:
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md),
[`sigmoid_length()`](https://sizespectrum.org/mizer/reference/sigmoid_length.md),
[`sigmoid_weight()`](https://sizespectrum.org/mizer/reference/sigmoid_weight.md)

## Examples

``` r
# Hump-shaped selectivity: rises from l25=10 to l50=15,
# then drops back to 50% at l50_right=40 and 25% at l25_right=50
sp <- list(a = 0.01, b = 3)
w <- c(1, 10, 100, 500, 1000)
double_sigmoid_length(w, l25 = 10, l50 = 15,
                      l50_right = 40, l25_right = 50,
                      species_params = sp)
#> [1] 0.09125627 0.24107143 0.71411741 0.58113240 0.33040411
```
