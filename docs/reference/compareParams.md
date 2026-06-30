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
#> Warning: For the species Sprat the value for `w_mat25` is not smaller than that of `w_mat`. I have corrected that by setting it to NA.
compareParams(params1, params2)
#> The following species parameters differ: Component “w_mat”: Mean absolute difference: 3, Component “w_mat25”: Mean absolute difference: 2.687875
#> 
#> The time_modified slots do not agree: Mean absolute difference: 111848.4
#> 
#> The maturity slots do not agree: Mean absolute difference: 0.04497162
#>   Max |diff|: Sprat: 0.568
#> 
#> The psi slots do not agree: Mean absolute difference: 0.06138336
#>   Max |diff|: Sprat: 0.403
#> 
```
