# Fishing effort used in simulation

Note that the array returned may not be exactly the same as the `effort`
argument that was passed in to
[`project()`](https://sizespectrum.org/mizer/reference/project.md). This
is because only the saved effort is stored (the frequency of saving is
determined by the argument `t_save`).

## Usage

``` r
getEffort(sim)
```

## Arguments

- sim:

  A MizerSim object

## Value

An array (time x gear) that contains the fishing effort by time and
gear.

## Examples

``` r
str(getEffort(NS_sim))
#>  num [1:44, 1:12] 0 0 0 0 0 0 0 0 0 0 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ time: chr [1:44] "1967" "1968" "1969" "1970" ...
#>   ..$ gear: chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
```
