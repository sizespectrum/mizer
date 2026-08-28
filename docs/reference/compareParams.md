# Compare two MizerParams objects and print out differences

Compare two MizerParams objects and print out differences

## Usage

``` r
compareParams(params1, params2, ...)
```

## Arguments

- params1:

  First MizerParams object

- params2:

  Second MizerParams object

- ...:

  Additional arguments passed to the method.

## Value

Invisibly returns a character vector of difference messages, one element
per difference. As a side effect, prints the differences in a
human-readable format.

## Examples

``` r
params1 <- NS_params
params2 <- params1
species_params(params2)$w_mat[1] <- 10
#> Warning: For the species Sprat the value for `w_mat25` is not smaller than that of `w_mat`. I have corrected that by marking it as missing so that its default will be used.
# Keep this example focused on the model parameter.
params2@time_modified <- params1@time_modified
compareParams(params1, params2)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> The following species parameters differ: Component "w_mat": Mean relative difference: 0.2307692, Component "w_mat25": Mean relative difference: 0.2307692
#> 
#> params2 has the following additional given species parameters: w_mat25
#> 
#> The following given species parameters differ: Component "w_mat": Mean relative difference: 0.2307692
#> 
#> The maturity slots do not agree: Mean absolute difference: 0.04497162
#>   Max |diff|: Sprat: 0.568
#> 
#> The psi slots do not agree: Mean absolute difference: 0.06138336
#>   Max |diff|: Sprat: 0.403
#> 
```
