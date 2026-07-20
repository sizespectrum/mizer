# S3 class for time x species arrays

Some functions in mizer return two-dimensional arrays (time x species)
holding quantities like biomass, abundance, or yield rate through time.
The `ArrayTimeBySpecies` class wraps these arrays to provide convenient
[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md), and
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md)
methods.

## Usage

``` r
ArrayTimeBySpecies(x, value_name = NULL, units = NULL, params = NULL)

is.ArrayTimeBySpecies(x)
```

## Arguments

- x:

  A matrix (time x species). For `is.ArrayTimeBySpecies()`, any object
  to test.

- value_name:

  A string giving the human-readable name for the value.

- units:

  A string giving the units (e.g. "g", "g/year").

- params:

  A `MizerParams` object holding the model that created the values.

## Value

An `ArrayTimeBySpecies` object (inherits from `matrix` and `array`).

`is.ArrayTimeBySpecies()` returns `TRUE` if `x` is an
`ArrayTimeBySpecies` object, `FALSE` otherwise.

## Details

An `ArrayTimeBySpecies` object behaves just like a regular matrix for
arithmetic operations and subsetting. It carries these lightweight
attributes:

- `value_name` – a human-readable name for the value (e.g. "Biomass").

- `units` – the units of the value (e.g. "g", "g/year").

- `params` – the `MizerParams` object that created the values.

## See also

[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md)

## Examples

``` r
# \donttest{
bio <- getBiomass(NS_sim)
is.ArrayTimeBySpecies(bio)
#> [1] TRUE
summary(bio)
#> Biomass [g] 
#> 44 times x 12 species
#> 
#>  Species          Min         Mean          Max
#>    Sprat 1.879468e+10 3.974130e+10 5.654305e+10
#>  Sandeel 7.342796e+11 2.099658e+12 3.737535e+12
#>   N.pout 2.528021e+11 3.254617e+11 3.929766e+11
#>  Herring 2.152333e+11 3.944814e+11 5.243289e+11
#>      Dab 1.200040e+10 1.649782e+10 1.953451e+10
#>  Whiting 1.983459e+11 2.210414e+11 2.605881e+11
#>     Sole 9.370624e+10 1.058612e+11 1.280000e+11
#>  Gurnard 1.067321e+11 1.291922e+11 1.457401e+11
#>   Plaice 1.950539e+12 2.260909e+12 2.756253e+12
#>  Haddock 5.358152e+11 6.424916e+11 8.205085e+11
#>      Cod 2.544990e+11 3.871450e+11 1.618389e+12
#>   Saithe 6.211381e+11 7.998952e+11 1.062987e+12
# }
```
