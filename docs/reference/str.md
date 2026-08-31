# Display the structure of mizer objects

Mizer provides `str()` methods for
[`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.md)
and [`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md)
objects, as well as
[`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md),
[`ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md)
and
[`ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)
objects. These methods produce a clean, compact overview of the object's
structure without polluting the console with large amounts of internal
data.

## Usage

``` r
# S3 method for class 'ArraySpeciesBySize'
str(object, ...)
# S3 method for class 'ArrayTimeBySpecies'
str(object, ...)
# S3 method for class 'ArrayTimeBySpeciesBySize'
str(object, ...)
# S3 method for class 'MizerSim'
str(object, max.level = NA, ...)
# S3 method for class 'MizerParams'
str(object, max.level = NA, ...)
```

## Arguments

- object:

  The object to display the structure of.

- max.level:

  Maximum level of nesting to print. Defaults to `NA` (no limit).

- ...:

  Further arguments. They are passed to the default `str()` method.

## Value

`NULL`, invisibly.

## See also

[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md),
[`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.md),
[`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md),
[`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md),
[`ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md),
[`ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)

## Examples

``` r
# \donttest{
str(NS_params)
#> Formal class 'MizerParams' [package "mizer"] with 0 slots
str(NS_sim)
#> Formal class 'MizerSim' [package "mizer"] with 0 slots
str(getEncounter(NS_params))
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#>  'ArraySpeciesBySize' num [1:12, 1:100] 0.299 0.453 0.502 0.575 0.492 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ sp: chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>  - attr(*, "value_name")= chr "Encounter rate"
#>  - attr(*, "units")= chr "g/year"
#>  - attr(*, "type")= chr "value"
#>  - attr(*, "representation")= chr "point"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
# }
```
