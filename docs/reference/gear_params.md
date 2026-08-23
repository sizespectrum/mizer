# Gear parameters

These functions allow you to get or set the gear parameters stored in a
MizerParams object. These are used by
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
to set up the selectivity and catchability and thus together with the
fishing effort determine the fishing mortality.

## Usage

``` r
gear_params(object)

gear_params(object) <- value

is.gear_params(x)
```

## Arguments

- object:

  A MizerParams object, a MizerSim object or a data frame

- value:

  A data frame with the new gear parameters.

- x:

  An object to test with `is.gear_params()`.

## Value

Data frame with gear parameters

`is.gear_params()` returns `TRUE` if `x` is a `gear_params` object,
`FALSE` otherwise.

## Details

The `gear_params` data has one row for each gear-species pair and one
column for each parameter that determines how that gear interacts with
that species. The columns are:

- `species` The name of the species

- `gear` The name of the gear

- `catchability` A number specifying how strongly this gear selects this
  species.

- `sel_func` The name of the function that calculates the selectivity
  curve.

- One column for each selectivity parameter needed by the selectivity
  functions.

For the details see
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md).

There can optionally also be a column `yield_observed` that allows you
to specify for each gear and species the total annual fisheries yield in
grams per year. This is used by
[`plotYieldObservedVsModel()`](https://sizespectrum.org/mizer/reference/plotYieldObservedVsModel.md),
which adds the yields up over the gears to get the observed yield of
each species, see
[`get_yield_observed()`](https://sizespectrum.org/mizer/reference/get_yield_observed.md).

The fishing effort, which is also needed to determine the fishing
mortality exerted by a gear is not set via the `gear_params` data frame
but is set with
[`initial_effort()`](https://sizespectrum.org/mizer/reference/initial_effort.md)
or is specified when calling
[`project()`](https://sizespectrum.org/mizer/reference/project.md).

If you change a gear parameter, this will be used to recalculate the
`selectivity` and `catchability` arrays by calling
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md),
unless you have previously set these by hand.

`gear_params<-` automatically sets the row names to contain the species
name and the gear name, separated by a comma and a space. The last
example below illustrates how this facilitates changing an individual
gear parameter.

## See also

Other functions for setting parameters:
[`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md),
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md),
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)

## Examples

``` r
params <- NS_params

# gears set up in example
gear_params(params)
#> An object of class "gear_params" containing 12 gear-species pairs for 4 gears:
#>        gear species   sel_func catchability
#>  Industrial   Sprat knife_edge            1
#>  Industrial Sandeel knife_edge            1
#>  Industrial  N.pout knife_edge            1
#>     Pelagic Herring knife_edge            1
#>        Beam     Dab knife_edge            1
#>       Otter Whiting knife_edge            1
#>        Beam    Sole knife_edge            1
#>       Otter Gurnard knife_edge            1
#>        Beam  Plaice knife_edge            1
#>       Otter Haddock knife_edge            1
#>       Otter     Cod knife_edge            1
#>       Otter  Saithe knife_edge            1
#> With 1 other parameters: knife_edge_size 

# setting totally different gears
gear_params(params) <- data.frame(
    gear = c("gear1", "gear2", "gear1"),
    species = c("Cod", "Cod", "Haddock"),
    catchability = c(0.5, 2, 1),
    sel_func = c("sigmoid_weight", "knife_edge", "sigmoid_weight"),
    sigmoidal_weight = c(1000, NA, 800),
    sigmoidal_sigma = c(100, NA, 100),
    knife_edge_size = c(NA, 1000, NA)
    )
gear_params(params)
#> An object of class "gear_params" containing 3 gear-species pairs for 2 gears:
#>   gear species       sel_func catchability
#>  gear1     Cod sigmoid_weight          0.5
#>  gear2     Cod     knife_edge          2.0
#>  gear1 Haddock sigmoid_weight          1.0
#> With 3 other parameters: sigmoidal_weight, sigmoidal_sigma, knife_edge_size 

# changing an individual entry
gear_params(params)["Cod, gear1", "catchability"] <- 0.8
```
