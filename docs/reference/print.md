# Print mizer objects

Mizer supplies `print()` methods for the array-like objects returned by
many rate and summary functions. These methods print a compact, readable
overview instead of the full matrix or array. The printed output reports
the value name, dimensions, units if known, and per-species minimum,
mean and maximum values.

## Arguments

- x:

  The object to print.

- ...:

  Further arguments. They are currently ignored by the mizer methods.

## Value

The printed object, invisibly.

## Details

For full numeric access, use the object itself as an ordinary matrix or
array or convert it to a long data frame with
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md).

## See also

[`summary()`](https://sizespectrum.org/mizer/reference/summary.md),
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md),
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md),
[`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md),
[`ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md),
[`ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)

## Examples

``` r
# \donttest{
enc <- getEncounter(NS_params)
print(enc)
#> Encounter rate (12 species x 100 sizes) [g/year] 
#>   Sprat: min=0.299 mean=2930 max=39600
#>   Sandeel: min=0.453 mean=3770 max=45500
#>   N.pout: min=0.502 mean=16800 max=148000
#>   Herring: min=0.575 mean=6240 max=80400
#>   Dab: min=0.492 mean=24000 max=267000
#>   Whiting: min=0.436 mean=17300 max=138000
#>   Sole: min=0.365 mean=12100 max=149000
#>   Gurnard: min=0.312 mean=11300 max=135000
#>   Plaice: min=0.232 mean=11800 max=113000
#>   Haddock: min=0.596 mean=24900 max=335000
#>   Cod: min=0.966 mean=52600 max=437000
#>   Saithe: min=0.771 mean=16400 max=188000

biomass <- getBiomass(NS_sim)
print(biomass)
#> Biomass (44 times x 12 species) [g] 
#>   Sprat: min=1.88e+10 mean=3.97e+10 max=5.65e+10
#>   Sandeel: min=7.34e+11 mean=2.1e+12 max=3.74e+12
#>   N.pout: min=2.53e+11 mean=3.25e+11 max=3.93e+11
#>   Herring: min=2.15e+11 mean=3.94e+11 max=5.24e+11
#>   Dab: min=1.2e+10 mean=1.65e+10 max=1.95e+10
#>   Whiting: min=1.98e+11 mean=2.21e+11 max=2.61e+11
#>   Sole: min=9.37e+10 mean=1.06e+11 max=1.28e+11
#>   Gurnard: min=1.07e+11 mean=1.29e+11 max=1.46e+11
#>   Plaice: min=1.95e+12 mean=2.26e+12 max=2.76e+12
#>   Haddock: min=5.36e+11 mean=6.42e+11 max=8.21e+11
#>   Cod: min=2.54e+11 mean=3.87e+11 max=1.62e+12
#>   Saithe: min=6.21e+11 mean=8e+11 max=1.06e+12
# }
```
