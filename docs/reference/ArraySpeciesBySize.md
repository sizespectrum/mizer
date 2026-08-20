# S3 class for species x size rate arrays

Many functions in mizer return two-dimensional arrays (species x size)
holding rates like encounter rate, feeding level, growth rate, mortality
etc. The `ArraySpeciesBySize` class wraps these arrays to provide
convenient
[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md), and
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md)
methods.

## Usage

``` r
ArraySpeciesBySize(
  x,
  value_name = NULL,
  units = NULL,
  type = NULL,
  params = NULL,
  representation = c("point", "average")
)

is.ArraySpeciesBySize(x)
```

## Arguments

- x:

  A matrix (species x size). For `is.ArraySpeciesBySize()`, any object
  to test.

- value_name:

  A string giving the human-readable name for the value.

- units:

  A string giving the units (e.g. "g/year", "1/year").

- type:

  The kind of quantity the values are, see
  [array_types](https://sizespectrum.org/mizer/reference/array_types.md):
  `"value"` (the default) for a rate or an amount, `"density"` for an
  amount per gram of body weight, `"proportion"` for a fraction. This is
  what tells
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) to
  multiply a density by the appropriate Jacobian when it is plotted
  against a length axis (`size_axis = "l"`), and to show a proportion
  against the whole of the interval from 0 to 1. The default, `NULL`,
  treats a `value_name` of `"Number density"` or units of `"1/g"` as a
  density, the way mizer recognised one before this attribute existed.

- params:

  A `MizerParams` object. Used for species colours, linetypes, and size
  ranges in the
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) method.

- representation:

  Either `"point"` (the default) for a quantity sampled at the grid
  nodes, or `"average"` for a finite-volume bin average. A bin-averaged
  quantity is drawn at the geometric bin centre rather than the left bin
  edge, but only when the model uses second-order bin-averaging
  (`second_order_w[["bin_average"]]`), so default plots are unchanged.

## Value

An `ArraySpeciesBySize` object (inherits from `matrix` and `array`).

`is.ArraySpeciesBySize()` returns `TRUE` if `x` is an
`ArraySpeciesBySize` object, `FALSE` otherwise.

## Details

An `ArraySpeciesBySize` object behaves just like a regular matrix for
arithmetic operations and subsetting. It carries a few lightweight
attributes:

- `value_name` – a human-readable name for the value (e.g. "Encounter
  rate").

- `units` – the units of the rate (e.g. "g/year").

- `type` – the kind of quantity the values are.

## See also

[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md)

## Examples

``` r
# \donttest{
enc <- getEncounter(NS_params)
is.ArraySpeciesBySize(enc)
#> [1] TRUE
summary(enc)
#> Encounter rate [g/year] 
#> 12 species x 100 sizes
#> 
#>  Species       Min      Mean       Max
#>    Sprat 0.2992076  2929.178  39573.31
#>  Sandeel 0.4528175  3768.983  45507.81
#>   N.pout 0.5019776 16840.828 147886.62
#>  Herring 0.5752333  6241.503  80375.40
#>      Dab 0.4916095 24004.843 266704.99
#>  Whiting 0.4362525 17348.364 137539.48
#>     Sole 0.3646753 12087.308 148784.12
#>  Gurnard 0.3122260 11327.057 135351.68
#>   Plaice 0.2323659 11800.440 113242.16
#>  Haddock 0.5964130 24932.923 334719.58
#>      Cod 0.9658343 52646.610 436916.91
#>   Saithe 0.7709631 16377.321 187775.81
# }
```
