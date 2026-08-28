# Plotting observed vs. model yields

**\[experimental\]** If yield observations are available for at least
some species via the `yield_observed` column, this function plots the
yield of each species in the model against the observed yields. When
called with a MizerSim object, the plot will use the model yields
predicted for the final time step in the simulation.

## Usage

``` r
plotYieldObservedVsModel(
  object,
  species = NULL,
  ratio = FALSE,
  log_scale = TRUE,
  return_data = FALSE,
  labels = TRUE,
  show_unobserved = FALSE,
  gear = NULL,
  ...
)
```

## Arguments

- object:

  An object of class
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  or
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md).

- species:

  The species to be included. Optional. By default all observed yields
  will be included. A vector of species names, or a numeric vector with
  the species indices, or a logical vector indicating for each species
  whether it is to be included (TRUE) or not.

- ratio:

  Whether to plot model yield vs. observed yield (FALSE) or the ratio of
  model : observed yield (TRUE). Default is FALSE.

- log_scale:

  Whether to plot on the log10 scale (TRUE) or not (FALSE). For the
  non-ratio plot this applies for both axes, for the ratio plot only the
  x-axis is on the log10 scale. Default is TRUE.

- return_data:

  Whether to return the data frame for the plot (TRUE) or not (FALSE).
  Default is FALSE.

- labels:

  Whether to show text labels for each species (TRUE) or not (FALSE).
  Default is TRUE.

- show_unobserved:

  Whether to include also species for which no yield observation is
  available. If TRUE, these species will be shown as if their observed
  yield was equal to the model yield.

- gear:

  The gears to be included. Optional. By default the catch of all gears
  is included. A vector of gear names. Only species caught by the
  selected gears are shown.

- ...:

  For `plotlyYieldObservedVsModel()`, additional arguments passed to
  [`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md).
  Otherwise unused.

## Value

A ggplot2 object with the plot of model yield by species compared to
observed yield. If `return_data = TRUE`, the data frame used to create
the plot is returned instead of the plot.

## Details

Before you can use this function you will need to have added a
`yield_observed` column to your model which gives the observed yield in
grams per year. Its home is the gear parameter data frame, see
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
where you give the yield for each gear-species pair and this function
adds them up over the gears. For backwards compatibility a
`yield_observed` column in the species parameter data frame is also
accepted, see
[`get_yield_observed()`](https://sizespectrum.org/mizer/reference/get_yield_observed.md).
For species for which you have no observed yield, you should set the
value in the `yield_observed` column to 0 or NA.

If a species is caught by several gears, both the model yield and the
observed yield are summed over the gears. With the `gear` argument you
can restrict the comparison to a subset of the gears, in which case only
the catch of those gears enters on both axes. Because the species
parameter data frame only holds the yield summed over all gears, the
observations then have to come from the gear parameters.

The total relative error is shown in the caption of the plot, calculated
by \$\$TRE = \sum_i\|1-\rm{ratio_i}\|\$\$ where \\\rm{ratio_i}\\ is the
ratio of model yield / observed yield for species i.

## Examples

``` r
# create an example
params <- NS_params
# In this model each species is caught by a single gear, so there is one
# row in the gear parameters for each species, in the same order.
# Species without an observation get NA.
gear_params(params)$yield_observed <-
    c(NA, NA, NA, 3e11, 4e9, 4e10, 5e10, NA, 2e11, 6e10, 3e11, NA)

# Plot with default options
plotYieldObservedVsModel(params)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> The following species are not being fished in your model and will not be included in the plot: Sprat, Sandeel, N.pout.


# Plot including also species without observations
plotYieldObservedVsModel(params, show_unobserved = TRUE)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> The following species are not being fished in your model and will not be included in the plot: Sprat, Sandeel, N.pout.


# Show the ratio instead
plotYieldObservedVsModel(params, ratio = TRUE)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> The following species are not being fished in your model and will not be included in the plot: Sprat, Sandeel, N.pout.


# If several gears catch the same species, their yields are added up.
# Give Cod a second gear that takes a quarter of the observed yield.
gp <- gear_params(params)
gp["Cod, Otter", "yield_observed"] <- 3e11 * 0.75
extra <- gp["Cod, Otter", ]
extra$gear <- "Gillnet"
extra$yield_observed <- 3e11 * 0.25
gear_params(params) <- rbind(gp, extra)

# Compare only the catch of the Otter gear against its observation
plotYieldObservedVsModel(params, gear = "Otter")
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> The following species are not being caught by the selected gear and will not be included in the plot: Sprat, Sandeel, N.pout, Herring, Dab, Sole, Plaice.
```
