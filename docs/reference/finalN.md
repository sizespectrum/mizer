# Size spectra at end of simulation

Size spectra at end of simulation

## Usage

``` r
finalN(sim)

finalNResource(sim)

idxFinalT(sim)
```

## Arguments

- sim:

  A MizerSim object

## Value

For `finalN()`: An `ArraySpeciesBySize` object (species x size) holding
the consumer number densities at the end of the simulation

For `finalNResource()`: A vector holding the resource number densities
at the end of the simulation for all size classes

For `idxFinalT()`: An integer giving the index for extracting the
results for the final time step

## Examples

``` r
str(finalN(NS_sim))
#>  'ArraySpeciesBySize' num [1:12, 1:100] 1.53e+13 5.89e+12 1.10e+14 1.37e+13 1.16e+11 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ sp: chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>  - attr(*, "value_name")= chr "Number density"
#>  - attr(*, "representation")= chr "point"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots

# This could also be obtained using `N()` and `idxFinalT()`
identical(N(NS_sim)[idxFinalT(NS_sim), , ], finalN(NS_sim))
#> [1] FALSE
str(finalNResource(NS_sim))
#>  'ArrayResourceBySize' Named num [1:218] 4.88e+35 3.40e+35 2.36e+35 1.64e+35 1.14e+35 ...
#>  - attr(*, "names")= chr [1:218] "8.73e-13" "1.04e-12" "1.24e-12" "1.48e-12" ...
#>  - attr(*, "value_name")= chr "Number density"
#>  - attr(*, "units")= chr "1/g"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
idx <- idxFinalT(NS_sim)
idx
#> [1] 44
# This coincides with
length(getTimes(NS_sim))
#> [1] 44
# and corresponds to the final time
getTimes(NS_sim)[idx]
#> [1] 2010
# We can use this index to extract the result at the final time
identical(N(NS_sim)[idx, , ], finalN(NS_sim))
#> [1] FALSE
identical(NResource(NS_sim)[idx, ], finalNResource(NS_sim))
#> [1] TRUE
```
