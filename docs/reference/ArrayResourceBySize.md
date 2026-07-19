# S3 class for resource size spectra

Several functions in mizer return a vector over the full size grid
holding a resource-related quantity such as the resource number density,
the resource mortality, the intrinsic resource birth rate or carrying
capacity. The `ArrayResourceBySize` class wraps these vectors to provide
convenient
[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md), and
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md)
methods.

## Usage

``` r
ArrayResourceBySize(x, value_name = NULL, units = NULL, params = NULL)

is.ArrayResourceBySize(x)
```

## Arguments

- x:

  A numeric vector over the full size grid. For
  `is.ArrayResourceBySize()`, any object to test.

- value_name:

  A string giving the human-readable name for the value.

- units:

  A string giving the units (e.g. "1/year").

- params:

  A `MizerParams` object. Used for the resource colour and the size grid
  in the [`plot()`](https://sizespectrum.org/mizer/reference/plot.md)
  method.

## Value

An `ArrayResourceBySize` object (inherits from `numeric`).

`is.ArrayResourceBySize()` returns `TRUE` if `x` is an
`ArrayResourceBySize` object, `FALSE` otherwise.

## Details

An `ArrayResourceBySize` object behaves just like a regular numeric
vector for arithmetic operations and subsetting. It carries three
lightweight attributes:

- `value_name` – a human-readable name for the value (e.g. "Resource
  mortality").

- `units` – the units of the value (e.g. "1/year").

- `params` – the `MizerParams` object that the value was computed from.

## See also

[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md)

## Examples

``` r
# \donttest{
mort <- getResourceMort(NS_params)
is.ArrayResourceBySize(mort)
#> [1] TRUE
summary(mort)
#> Resource mortality [1/year] 
#> 226 sizes
#> 
#>           Min     Mean      Max
#>  1.228364e-16 6.936203 24.10931
plot(mort)

# }
```
