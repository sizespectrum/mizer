# Gear parameters

These functions allow you to get or set the gear parameters stored in a
MizerParams object. These are used by
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
to set up the selectivity and catchability and thus together with the
fishing effort determine the fishing mortality.

## Usage

``` r
gear_params(params)

gear_params(params) <- value
```

## Arguments

- params:

  A MizerParams object

- value:

  A data frame with the gear parameters.

## Value

Data frame with gear parameters

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
to specify for each gear and species the total annual fisheries yield.

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

[`validGearParams()`](https://sizespectrum.org/mizer/reference/validGearParams.md)

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
#>                           gear species   sel_func knife_edge_size catchability
#> Sprat, Industrial   Industrial   Sprat knife_edge              13            1
#> Sandeel, Industrial Industrial Sandeel knife_edge               4            1
#> N.pout, Industrial  Industrial  N.pout knife_edge              23            1
#> Herring, Pelagic       Pelagic Herring knife_edge              99            1
#> Dab, Beam                 Beam     Dab knife_edge              21            1
#> Whiting, Otter           Otter Whiting knife_edge              75            1
#> Sole, Beam                Beam    Sole knife_edge              78            1
#> Gurnard, Otter           Otter Gurnard knife_edge              39            1
#> Plaice, Beam              Beam  Plaice knife_edge             105            1
#> Haddock, Otter           Otter Haddock knife_edge             165            1
#> Cod, Otter               Otter     Cod knife_edge            1606            1
#> Saithe, Otter            Otter  Saithe knife_edge            1076            1

# setting totally different gears
gear_params(params) <- data.frame(
    gear = c("gear1", "gear2", "gear1"),
    species = c("Cod", "Cod", "Haddock"),
    catchability = c(0.5, 2, 1),
    sel_fun = c("sigmoid_weight", "knife_edge", "sigmoid_weight"),
    sigmoidal_weight = c(1000, NA, 800),
    sigmoidal_sigma = c(100, NA, 100),
    knife_edge_size = c(NA, 1000, NA)
    )
gear_params(params)
#>                 gear species catchability        sel_fun sigmoidal_weight
#> Cod, gear1     gear1     Cod          0.5 sigmoid_weight             1000
#> Cod, gear2     gear2     Cod          2.0     knife_edge               NA
#> Haddock, gear1 gear1 Haddock          1.0 sigmoid_weight              800
#>                sigmoidal_sigma knife_edge_size   sel_func
#> Cod, gear1                 100            1606 knife_edge
#> Cod, gear2                  NA            1000 knife_edge
#> Haddock, gear1             100             165 knife_edge

# changing an individual entry
gear_params(params)["Cod, gear1", "catchability"] <- 0.8
```
