# S3 class for time x species x size arrays

Some functions in mizer return three-dimensional arrays (time x species
x size) holding quantities like fishing mortality, feeding level, or
predation mortality through time. The `ArrayTimeBySpeciesBySize` class
wraps these arrays to provide convenient
[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md),
[`animate()`](https://sizespectrum.org/mizer/reference/animate.md), and
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md)
methods.

## Usage

``` r
ArrayTimeBySpeciesBySize(
  x,
  value_name = NULL,
  units = NULL,
  type = NULL,
  params = NULL,
  representation = c("point", "average")
)

is.ArrayTimeBySpeciesBySize(x)
```

## Arguments

- x:

  A 3D array (time x species x size). For
  `is.ArrayTimeBySpeciesBySize()`, any object to test.

- value_name:

  A string giving the human-readable name for the value.

- units:

  A string giving the units (e.g. "1/year").

- type:

  The kind of quantity the values are, see
  [`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md)
  and
  [array_types](https://sizespectrum.org/mizer/reference/array_types.md).

- params:

  A `MizerParams` object. Used for species colours, linetypes, and size
  ranges in the
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) and
  [`animateSpectra()`](https://sizespectrum.org/mizer/reference/animate.md)
  methods.

- representation:

  Either `"point"` (the default) for a quantity sampled at the grid
  nodes, or `"average"` for a finite-volume bin average, which is then
  drawn at the geometric bin centre when the model uses second-order
  bin-averaging (`second_order_w[["bin_average"]]`). See
  [`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md).

## Value

An `ArrayTimeBySpeciesBySize` object (inherits from `array`).

`is.ArrayTimeBySpeciesBySize()` returns `TRUE` if `x` is an
`ArrayTimeBySpeciesBySize` object, `FALSE` otherwise.

## Details

An `ArrayTimeBySpeciesBySize` object behaves just like a regular array
for arithmetic operations and subsetting. It carries these lightweight
attributes:

- `value_name` – a human-readable name for the value (e.g. "Fishing
  mortality").

- `units` – the units of the value (e.g. "1/year").

- `type` – the kind of quantity the values are.

- `params` – the `MizerParams` object that the value was computed from.

## See also

[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md),
[`animateSpectra()`](https://sizespectrum.org/mizer/reference/animate.md)

## Examples

``` r
# \donttest{
fmort <- getFMort(NS_sim)
is.ArrayTimeBySpeciesBySize(fmort)
#> [1] TRUE
summary(fmort)
#> Fishing mortality [1/year] 
#> 44 times x 12 species x 100 sizes
#> 
#>  Species          Min       Mean       Max
#>    Sprat 0.0000000000 0.12879552 2.1827923
#>  Sandeel 0.0000000000 0.11611087 1.3076336
#>   N.pout 0.0000000000 0.13718398 2.1228994
#>  Herring 0.0077894362 0.18281500 1.4419293
#>      Dab 0.0017746971 0.02322755 0.1633433
#>  Whiting 0.0126246154 0.17579846 1.3138868
#>     Sole 0.0268320133 0.15939358 1.0259979
#>  Gurnard 0.0000000000 0.00541050 0.1054407
#>   Plaice 0.0089964121 0.18786555 0.8670560
#>  Haddock 0.0015854377 0.22672527 1.4275500
#>      Cod 0.0494762790 0.36512704 1.0721072
#>   Saithe 0.0009633361 0.15208058 1.2032857
plot(fmort, time = 2007)

# }
```
