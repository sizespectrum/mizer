# S3 class for time x resource-size arrays

The [`NResource()`](https://sizespectrum.org/mizer/reference/N.md)
function returns a two-dimensional array (time x size) holding the
resource number density through time. The `ArrayTimeByResourceBySize`
class wraps this array to provide convenient
[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md), and
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md)
methods.

## Usage

``` r
ArrayTimeByResourceBySize(x, value_name = NULL, units = NULL, params = NULL)
```

## Arguments

- x:

  A matrix (time x size).

- value_name:

  A string giving the human-readable name for the value.

- units:

  A string giving the units (e.g. "1/g").

- params:

  A `MizerParams` object. Used for the resource colour and the size grid
  in the [`plot()`](https://sizespectrum.org/mizer/reference/plot.md)
  method.

## Value

An `ArrayTimeByResourceBySize` object (inherits from `matrix` and
`array`).

## Details

An `ArrayTimeByResourceBySize` object behaves just like a regular matrix
for arithmetic operations and subsetting. It carries these lightweight
attributes:

- `value_name` – a human-readable name for the value (e.g. "Number
  density").

- `units` – the units of the value (e.g. "1/g").

- `params` – the `MizerParams` object that the value was computed from.

## See also

[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md)

## Examples

``` r
# \donttest{
nr <- NResource(NS_sim)
is.ArrayTimeByResourceBySize(nr)
#> [1] TRUE
summary(nr)
#> Number density [1/g] 
#> 44 times x 218 sizes
#> 
#>  Min         Mean          Max
#>    0 7.361221e+33 4.878226e+35
plot(nr)
#> Warning: log-10 transformation introduced infinite values.

# }
```
