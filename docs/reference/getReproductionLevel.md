# Get reproduction level

The reproduction level is the ratio between the density-dependent
reproduction rate and the maximal reproduction rate.

## Usage

``` r
getReproductionLevel(params)
```

## Arguments

- params:

  A MizerParams object

## Value

A named vector with the reproduction level for each species.

## Examples

``` r
getReproductionLevel(NS_params)
#>      Sprat    Sandeel     N.pout    Herring        Dab    Whiting       Sole 
#> 0.99074238 0.99987053 0.92829319 0.99198802 0.99578514 0.98718674 0.99643774 
#>    Gurnard     Plaice    Haddock        Cod     Saithe 
#> 0.44189813 0.08022106 0.94443443 0.99993658 0.99767830 

# The reproduction level can be changed without changing the steady state:
params <- setBevertonHolt(NS_params, reproduction_level = 0.9)
#> Warning: The following species require an unrealistic value greater than 1 for `erepro`: Gurnard, Plaice
getReproductionLevel(params)
#>   Sprat Sandeel  N.pout Herring     Dab Whiting    Sole Gurnard  Plaice Haddock 
#>     0.9     0.9     0.9     0.9     0.9     0.9     0.9     0.9     0.9     0.9 
#>     Cod  Saithe 
#>     0.9     0.9 

# The result is the ratio of RDD and R_max
identical(getRDD(params) / species_params(params)$R_max,
          getReproductionLevel(params))
#> [1] TRUE
```
