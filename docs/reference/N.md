# Time series of size spectra

Fetch the simulation results for the size spectra over time.

## Usage

``` r
N(sim)

NResource(sim)
```

## Arguments

- sim:

  A MizerSim object

## Value

For `N()`: An `ArrayTimeBySpeciesBySize` object (time x species x size)
with the number density of consumers.

For `NResource()`: An array (time x size) with the number density of
resource

## Examples

``` r
str(N(NS_sim))
#>  'ArrayTimeBySpeciesBySize' num [1:44, 1:12, 1:100] 1.67e+13 1.67e+13 1.69e+13 1.66e+13 1.65e+13 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ time: chr [1:44] "1967" "1968" "1969" "1970" ...
#>   ..$ sp  : chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w   : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>  - attr(*, "value_name")= chr "Number density"
#>  - attr(*, "units")= chr "1/g"
#>  - attr(*, "representation")= chr "point"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
str(NResource(NS_sim))
#>  'ArrayTimeByResourceBySize' num [1:44, 1:218] 4.88e+35 4.88e+35 4.88e+35 4.88e+35 4.88e+35 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ time: chr [1:44] "1967" "1968" "1969" "1970" ...
#>   ..$ w   : chr [1:218] "8.73e-13" "1.04e-12" "1.24e-12" "1.48e-12" ...
#>  - attr(*, "value_name")= chr "Number density"
#>  - attr(*, "units")= chr "1/g"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
```
