# Calculate selectivity from gear parameters

This function calculates the selectivity for each gear, species and size
from the gear parameters. It is called by
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
when the `selectivity` is not set by the user. The returned array is
initialised to zero, so gear-species combinations that are not listed in
`gear_params(params)` remain zero. For each listed combination the
function named in `sel_func` is called with `w = params@w`, the
corresponding species parameters, and the selectivity parameters from
the matching row in `gear_params(params)`.

## Usage

``` r
calc_selectivity(params)
```

## Arguments

- params:

  A MizerParams object

## Value

An array (gear x species x size) with the selectivity values

## Bin-averaged selectivity

By default the selectivity is point-sampled at the grid nodes
`params@w`, i.e. at the left edge of each size bin. This is only
first-order accurate in the bin size when the selectivity is used in the
finite-volume update of the size spectrum. When the `bin_average` entry
of the
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
slot is `TRUE`, each selectivity function is instead integrated over its
size bin, so that `selectivity[g, i, j]` holds the bin average \$\$\bar
S\_{g,i,j} = \frac{1}{\Delta w_j} \int\_{w_j}^{w\_{j+1}} S\_{g,i}(w)\\
dw.\$\$ The integral is evaluated with a composite-midpoint rule on a
log-spaced sub-grid of each bin, mirroring the bin-integrated predation
kernel. This lifts the fishing mortality towards second order at no
extra runtime cost (the integration happens once here, the rate
functions are unchanged). A welcome side effect is that a knife-edge
gear then gets the exact fraction of the straddling bin that lies above
the knife edge, removing a grid artefact.

## Examples

``` r
params <- NS_params
str(calc_selectivity(params))
#>  num [1:4, 1:12, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ gear: chr [1:4] "Industrial" "Pelagic" "Beam" "Otter"
#>   ..$ sp  : chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w   : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
calc_selectivity(params)["Pelagic", "Herring", ]
#>   0.001 0.00119 0.00142  0.0017 0.00203 0.00242 0.00289 0.00345 0.00411 0.00491 
#>       0       0       0       0       0       0       0       0       0       0 
#> 0.00586 0.00699 0.00834 0.00995  0.0119  0.0142  0.0169  0.0202  0.0241  0.0288 
#>       0       0       0       0       0       0       0       0       0       0 
#>  0.0343  0.0409  0.0489  0.0583  0.0696   0.083  0.0991   0.118   0.141   0.168 
#>       0       0       0       0       0       0       0       0       0       0 
#>   0.201    0.24   0.286   0.342   0.408   0.486    0.58   0.693   0.827   0.987 
#>       0       0       0       0       0       0       0       0       0       0 
#>    1.18     1.4    1.68       2    2.39    2.85     3.4    4.06    4.84    5.78 
#>       0       0       0       0       0       0       0       0       0       0 
#>     6.9    8.23    9.82    11.7      14    16.7    19.9    23.8    28.4    33.8 
#>       0       0       0       0       0       0       0       0       0       0 
#>    40.4    48.2    57.5    68.7    81.9    97.8     117     139     166     198 
#>       0       0       0       0       0       0       1       1       1       1 
#>     237     282     337     402     480     573     683     816     973    1160 
#>       1       1       1       1       1       1       1       1       1       1 
#>    1390    1650    1970    2360    2810    3350    4000    4780    5700    6800 
#>       1       1       1       1       1       1       1       1       1       1 
#>    8120    9690   11600   13800   16500   19600   23400   28000   33400   39900 
#>       1       1       1       1       1       1       1       1       1       1 
```
