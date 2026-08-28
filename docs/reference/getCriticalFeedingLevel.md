# Get critical feeding level

The critical feeding level is the feeding level at which the food intake
is just high enough to cover the metabolic costs, with nothing left over
for growth or reproduction.

## Usage

``` r
getCriticalFeedingLevel(params)
```

## Arguments

- params:

  A MizerParams object

## Value

An `ArraySpeciesBySize` object (species x size) with the critical
feeding level

## Examples

``` r
# \donttest{
str(getFeedingLevel(NS_params))
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#>  'ArraySpeciesBySize' num [1:12, 1:100] 0.622 0.622 0.605 0.621 0.616 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ sp: chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
#>  - attr(*, "value_name")= chr "Feeding level"
#>  - attr(*, "type")= chr "proportion"
#>  - attr(*, "representation")= chr "point"
#>  - attr(*, "params")=Formal class 'MizerParams' [package "mizer"] with 48 slots
# }
```
