# Convert mizer arrays to data frames

The `as.data.frame()` methods for mizer array classes turn matrix- and
array-like results into tidy long-form data frames, with one row per
observed combination of species, size and/or time. The numeric result is
always stored in a column called `value`.

## Usage

``` r
# S3 method for class 'ArraySpeciesBySize'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
# S3 method for class 'ArrayTimeBySpecies'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
# S3 method for class 'ArrayTimeBySpeciesBySize'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `ArraySpeciesBySize`, `ArrayTimeBySpecies` or
  `ArrayTimeBySpeciesBySize` object.

- row.names:

  Optional and included only for compatibility with the base generic.
  `NULL` or a character vector giving the row names for the data frame.

- optional:

  Optional and included only for compatibility with the base generic. A
  logical value. If `TRUE`, setting row names and converting column
  names (to syntactic names) is optional.

- ...:

  Further arguments. They are currently ignored by the mizer methods.

## Value

A data frame in long format.

## Details

The returned columns are:

- `ArraySpeciesBySize`: `w`, `value`, `Species`.

- `ArrayTimeBySpecies`: `time`, `value`, `Species`.

- `ArrayTimeBySpeciesBySize`: `time`, `Species`, `w`, `value`.

If the original object has non-numeric or missing dimension names,
sequential indices are used for the `time` or `w` columns. Species names
are taken from the row, column or dimension names of the original
object.

## See also

[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`str()`](https://sizespectrum.org/mizer/reference/str.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md),
[`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md),
[`ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md),
[`ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)

## Examples

``` r
# \donttest{
enc <- getEncounter(NS_params)
head(as.data.frame(enc))
#>       w     value Species
#> 1 0.001 0.2992076   Sprat
#> 2 0.001 0.4528175 Sandeel
#> 3 0.001 0.5019776  N.pout
#> 4 0.001 0.5752333 Herring
#> 5 0.001 0.4916095     Dab
#> 6 0.001 0.4362525 Whiting

biomass <- getBiomass(NS_sim)
head(as.data.frame(biomass))
#>   time       value Species
#> 1 1967 50836384568   Sprat
#> 2 1968 55698644563   Sprat
#> 3 1969 54815404811   Sprat
#> 4 1970 53170193699   Sprat
#> 5 1971 51643049443   Sprat
#> 6 1972 49812211590   Sprat
# }
```
