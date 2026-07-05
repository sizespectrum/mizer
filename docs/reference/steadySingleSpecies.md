# Set initial abundances to solution of steady-state equation with current rates

**\[experimental\]** This first calculates growth and death rates that
arise from the current initial abundances. Then it solves the
steady-state equation with these growth and death rates and the current
abundance at the smallest size. It sets the initial abundances of the
selected species to this solution.

## Usage

``` r
steadySingleSpecies(
  params,
  species = NULL,
  keep = c("egg", "biomass", "number")
)
```

## Arguments

- params:

  A MizerParams object

- species:

  The species to be selected. Optional. By default all target species
  are selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- keep:

  A string determining which quantity is to be kept constant. The
  choices are "egg" which keeps the egg density constant, "biomass"
  which keeps the total biomass of the species constant and "number"
  which keeps the total number of individuals constant.

## Value

A MizerParams object in which the initial abundances of the selected
species are changed to their single-species steady state abundances.

## Details

The function only changes the initial abundances. It does not adjust the
reproduction parameters or any other parameters. Therefore the result of
applying this function is of course not a steady state, because after
changing the abundances of the selected species the growth, death and
reproduction rates will have changed.

If the `keep` argument is supplied, the solution for the selected
species are rescaled to keep the specified quantity at the value they
had before calling this function.

## Examples

``` r
# Set initial abundance of Cod to its single-species steady state
params <- steadySingleSpecies(NS_params, species = "Cod")
```
