# The Trait-Based Model

## Introduction

The trait-based size spectrum model can be derived as a simplification
of the general model outline in the [model
description](https://sizespectrum.org/mizer/articles/model_description.md)
section. It is more complicated than a [community
model](https://sizespectrum.org/mizer/articles/community_model.md) and
the most significant difference between the two is that while the
community model aggregates all species into a single spectrum, the
trait-based model resolves many species.

In a trait-based model the maximum size is considered to be the most
important *trait* characterizing a species. Most of the species-specific
parameters, such as \\\beta\\ and \\\sigma\\, are the same for all
species. Other model parameters are determined by the maximum size. For
example, the weight at maturation is a fixed fraction of the maximum
size. The maximum sizes of the species are spread evenly on a
logarithmic scale. The number of species is not important and does not
affect the general dynamics of the model.

## Setting up a trait-based model

To help set up a trait-based model, there is a wrapper function,
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md).
Like the
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
function described in the section on the [community
model](https://sizespectrum.org/mizer/articles/community_model.md), this
function can take many arguments. Most of them have default values so
you don’t need to worry about them for the moment. See the help page for
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)
for more details.

One of the key differences between the community type model and the
trait-based model is that reproduction and egg production are
considered. In the community model, reproduction is constant and there
is no relationship between the abundance in the community and egg
production. In the trait-based model, the egg production is dependent on
mature individuals investing part of their energy income into
reproduction. The relationship between the energy invested into
reproduction and the actual rate of egg production is modelled using a
*Beverton-Holt* type function (the default in `mizer`, see [the section
on density-dependence in
reproduction](https://sizespectrum.org/mizer/articles/model_description.html#density-dependence-in-reproduction))
where the reproduction rate \\R_i\\ (numbers per time) approaches a
maximum as the energy invested increases.

Here we set up the model to have 10 species, with maximum sizes ranging
from 10 g to 100 kg. All the other parameters have default values.

``` r
params <- newTraitParams(no_sp = 10, min_w_max = 10, max_w_max = 1e5)
```

This function returns an object of type `MizerParams` that holds all the
model information, including species parameters. This object can
therefore be interrogated in the same way as described in [the section
on the community
model](https://sizespectrum.org/mizer/articles/community_model.md).

``` r
summary(params)
```

    ## An object of class "MizerParams" 
    ## mizer version: 3.3.1
    ## Created: 2026-08-28 14:51:43
    ## Modified: 2026-08-28 14:51:44
    ## Consumer size spectrum:
    ##  minimum size:   0.001
    ##  maximum size:   1e+05
    ##  no. size bins:  161
    ## Resource size spectrum:
    ##  minimum size:   1e-10
    ##  maximum size:   2.23872
    ##  no. size bins:  208 (301 size bins in total)
    ## Steady state:
    ##  biomass drift:  0.25 /year  (not at steady state, largest in 1 - run tuneSteadyState())
    ## Species details:
    ## An object of class "species_params" containing parameters for 10 species:
    ##  species        w_inf        w_mat w_min  f0   fc beta
    ##        1 8.912509e+00     2.511886 0.001 0.6 0.25  100
    ##        2 2.511886e+01     7.079458 0.001 0.6 0.25  100
    ##        3 7.079458e+01    19.952623 0.001 0.6 0.25  100
    ##        4 1.995262e+02    56.234133 0.001 0.6 0.25  100
    ##        5 5.623413e+02   158.489319 0.001 0.6 0.25  100
    ##        6 1.584893e+03   446.683592 0.001 0.6 0.25  100
    ##        7 4.466836e+03  1258.925412 0.001 0.6 0.25  100
    ##        8 1.258925e+04  3548.133892 0.001 0.6 0.25  100
    ##        9 3.548134e+04 10000.000000 0.001 0.6 0.25  100
    ##       10 1.000000e+05 28183.829313 0.001 0.6 0.25  100
    ## With 1 other parameters: sigma 
    ## 
    ## Fishing gear details:
    ## Gear          Effort  Target species 
    ##  ----------------------------------
    ## knife_edge_gear 0.00   1, 2, 3, 4, 5, 6, 7, 8, 9, 10

The summary shows us that now we have 10 species in the model, with
maximum sizes ranging from \\8.9125094\\ to \\10^{5}\\. The rather
strange-looking values for the sizes is due to the fact that the size
classes are equally spaced on a logarithmic scale.

The summary lists those sizes under `w_inf`, the von Bertalanffy
asymptotic size. In a general mizer model that is not the same thing as
`w_max`, which is the upper boundary of the size grid and by default
sits 50% above `w_inf` to leave room for individuals that grow past the
asymptotic size. The trait-based model switches diffusion off, so growth
is deterministic and no individual ever grows beyond `w_repro_max`, the
size at which all available income is invested into reproduction. There
is then nothing to leave room for, and
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)
puts all three at the same value:

``` r
species_params(params)[, c("w_inf", "w_repro_max", "w_max")]
```

    ## An object of class "species_params" containing parameters for 10 species:
    ##         w_inf  w_repro_max        w_max
    ##  8.912509e+00 8.912509e+00 8.912509e+00
    ##  2.511886e+01 2.511886e+01 2.511886e+01
    ##  7.079458e+01 7.079458e+01 7.079458e+01
    ##  1.995262e+02 1.995262e+02 1.995262e+02
    ##  5.623413e+02 5.623413e+02 5.623413e+02
    ##  1.584893e+03 1.584893e+03 1.584893e+03
    ##  4.466836e+03 4.466836e+03 4.466836e+03
    ##  1.258925e+04 1.258925e+04 1.258925e+04
    ##  3.548134e+04 3.548134e+04 3.548134e+04
    ##  1.000000e+05 1.000000e+05 1.000000e+05

So in this article the “maximum size” of a species is unambiguous: it is
the asymptotic size, the size at which growth stops and the size grid
ends.

The size at maturity (`w_mat`) is linearly related to the maximum size.
Each species has the same preferred predator-prey mass ratio parameter
values (`beta` and `sigma`, see [the section on predator/prey mass
ratio](https://sizespectrum.org/mizer/articles/model_description.html#sec:ppmr)).
There are \\161\\ size bins in the community and \\301\\ size bins
including the resource spectrum. Ignore the summary section on fishing
gear for the moment. This is explained later.

## Running the trait-based model

As with the [community
model](https://sizespectrum.org/mizer/articles/community_model.md), we
can project the trait-based model through time using the
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
function. Here we project the model for 75 years without any fishing
(the `effort` argument is set to 0). We use the default initial
population abundances so there is no need to pass in any initial
population values (see [the guide on running a
simulation](https://sizespectrum.org/mizer/articles/guide-run-simulation.md)).

``` r
sim <- project(params, t_max = 75, effort = 0)
```

This results in a `MizerSim` object that contains the abundances of the
community and resource spectra through time, as well as the original
`MizerParams` object. As with the community model, we can get a quick
overview of the results of the simulation by calling the
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) method:

``` r
plot(sim)
```

![Summary plot of the trait-based model without
fishing.](trait_model_files/figure-html/plot_comm_sim_2-1.png)

The summary plot has the same panels as the one generated by the
community model, but here you can see that all the species in the
community are plotted. The panels show the situation in the final time
step of the simulation, apart from the biomass through time plot. As
this is a trait-based model where all species fully interact with each
other, the predation mortality and feeding level by size is the same for
each species. In this simulation we turned fishing off so the fishing
mortality is 0. The size-spectra show the abundances at size to be
evenly spaced by log of maximum size.

## Example of a trophic cascade with the trait-based model

As with the community model, it is possible to use the trait-based model
to simulate a trophic cascade. Again, we perform two simulations, one
with fishing and one without. We therefore need to consider how fishing
gears and selectivity have been set up by the
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)
function.

The default fishing selectivity function is a knife-edge function, which
only selects individuals larger than 1000 g. There is also only one
fishing gear in operation, and this selects all of the species. You can
see this if you call the
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md)
method on the `params` argument we set up above. At the bottom of the
summary there is a section on *Fishing gear details*. You can see that
there is only one gear, called `knife_edge_gear` and that it selects
species 1 to 10. To control the size at which individuals are selected
there is a `knife_edge_size` argument to the
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)
function. This has a default value of 1000 g.

In `mizer` it is possible to include more than one fishing gear in the
model and for different species to be caught by different gears. We will
ignore this for now, but will explore it further below when we introduce
an industrial fishery to the trait-based model.

To set up the trait-based model to have fishing we set up the
`MizerParams` object in exactly the same way as we did before but here
the `knife_edge_size` argument is explicitly passed in for clarity:

``` r
params_knife <- newTraitParams(no_sp = 10, min_w_max = 10, max_w_max = 1e5,
    knife_edge_size = 1000)
```

First we perform a simulation without fishing in the same way we did
above by setting the `effort` argument to 0:

``` r
sim0 <- project(params_knife, effort = 0, t_max = 75)
```

Now we simulate with fishing. Here, we use an effort of 0.75. As
mentioned in [the section on trophic cascades in the community
model](https://sizespectrum.org/mizer/articles/community_model.html#sec:trophic_cascade_comm_model),
the fishing mortality on a species is calculated as the product of
effort, catchability and selectivity (see [the guide on setting up
fishing](https://sizespectrum.org/mizer/articles/guide-set-up-fishing.md)
for more details). Selectivity ranges between 0 (not selected) and 1
(fully selected). The default value of catchability is 1. Therefore, in
this simulation the fishing mortality of a fully selected individual is
simply equal to the effort. This effort is constant throughout the
duration of the simulation (however, mizer does allow variable effort).

``` r
sim1 <- project(params_knife, effort = 0.75, t_max = 75)
```

Again, we can plot the summary of the fished community using the default
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) function.
The knife-edge selectivity at 1000 g can be clearly seen in the fishing
mortality panel:

``` r
plot(sim1)
```

![Summary plot of the trait-based model with
fishing.](trait_model_files/figure-html/plot_trait_fmort-1.png)

The biomass panel of that plot deserves a second look: an effort of 0.75
is severe enough that six of the ten species are driven to extinction.
Only species 1, 4, 5 and 6 survive. Species 4 and 5 in fact end up more
abundant than they were before fishing started, because they are too
small to be caught themselves and have lost the larger species that used
to eat them. The trophic cascade we look at next therefore plays out in
a community that fishing has also reshaped.

The trophic cascade can be explored by comparing the abundance densities
of all species at size when the community is fished and unfished. As
with the community model, we are interested in the abundances in the
final time step. We do not need to plot that ratio by hand, because
[`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md)
makes exactly this comparison for us. It plots the difference between
two spectra relative to their average, \\2(N_1(w) - N_0(w)) / (N_1(w) +
N_0(w))\\, which stays between \\-2\\ and \\2\\ however extreme the
change is. We ask for the community `total` and `highlight` it so that
it stands out from the individual species:

``` r
plotSpectraRelative(sim0, sim1, total = TRUE, resource = FALSE,
                    highlight = "Total")
```

![The relative difference between the fished and the unfished
trait-based model at each size. The thick black total line rises above
zero around 0.1 g, dips below it between 1 and 100 g and collapses above
1000 g. Several species lie flat at -2 because fishing has driven them
extinct.](trait_model_files/figure-html/plot_relative_comm_abund2-1.png)

The impact of fishing on individuals larger than 1000 g can be clearly
seen: the total drops to \\-2\\ there, meaning that those sizes have
essentially been emptied. This relieves the predation pressure on their
smaller prey (the preferred predator-prey size ratio is given by the
\\\beta\\ parameter, which is set to 100 by default), leading to an
increase in their abundance. This in turn increases the predation
mortality on *their* smaller prey, which reduces their abundance and so
on, giving the total line its alternating rises and dips as we move down
the spectrum. The lines lying flat at \\-2\\ belong to species that have
been driven to extinction. Not all six of them appear: two have declined
so far that their densities have underflowed to zero, and a zero cannot
be shown on the logarithmic scale the spectra are prepared on.

This impact can also be seen by looking at the predation mortality by
size. The predation mortalities are retrieved using the
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md)
function, applied to the state of each simulation at its final time
step.
[`finalParams()`](https://sizespectrum.org/mizer/reference/getParams.md)
gives us that state as a `MizerParams` object, carrying the resource
spectrum along with the fish abundances:

``` r
m2_no_fishing <- getPredMort(finalParams(sim0))
m2_with_fishing <- getPredMort(finalParams(sim1))
```

In the trait-based model every species experiences the same predation
mortality at a given size, so it is enough to look at one of them. We
ask for the first species and for `all.sizes`, so that the curve is
drawn over the whole size range rather than only over the sizes that
species 1 itself reaches. The two mortalities are compared with
[`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md), which
draws two compatible mizer arrays in one plot:

``` r
plot2(m2_no_fishing, m2_with_fishing, name1 = "Unfished", name2 = "Fished",
      species = 1, all.sizes = TRUE)
```

![Predation mortality against size in the unfished (solid) and fished
(dashed) trait-based model. Fishing lowers the mortality on the smallest
individuals and raises it by more than a factor of two around 5
g.](trait_model_files/figure-html/plot_relative_trait_m2-1.png)

The peak of the predation mortality has moved: the small individuals
that used to be eaten most heavily are now much safer, while individuals
around 5 g face more than twice the mortality they did before. That
shift is the mechanism behind the cascade in the abundances above.

## Setting up an industrial fishing gear

In this section we want to operate an *industrial* fishery. Industrial
fishing targets the small zooplanktivorous species that are typically
used for fishmeal production.

In the previous simulations we had only one fishing gear and it targeted
all the species in the community. This gear had a knife-edge selectivity
that only selected species larger than 1 kg. We can see that when we
look at the gear parameters

``` r
gear_params(params)
```

    ## An object of class "gear_params" containing 10 gear-species pairs for 1 gear:
    ##             gear species   sel_func catchability
    ##  knife_edge_gear       1 knife_edge            1
    ##  knife_edge_gear       2 knife_edge            1
    ##  knife_edge_gear       3 knife_edge            1
    ##  knife_edge_gear       4 knife_edge            1
    ##  knife_edge_gear       5 knife_edge            1
    ##  knife_edge_gear       6 knife_edge            1
    ##  knife_edge_gear       7 knife_edge            1
    ##  knife_edge_gear       8 knife_edge            1
    ##  knife_edge_gear       9 knife_edge            1
    ##  knife_edge_gear      10 knife_edge            1
    ## With 1 other parameters: knife_edge_size

We will expand the model to include multiple fishing gears. This
requires us to look more closely at how fishing gears are handled in
`mizer`. In `mizer` it is possible for a fishing gear to catch only a
subset of the species in the model. This is useful because when running
a simulation with
[`project()`](https://sizespectrum.org/mizer/reference/project.md) you
can specify the effort per gear and so you can turn gears on or off as
you want. Each gear has a selectivity curve for each species.

We will set up the model to include two fishing gears: an `industrial`
gear that only catches species with a maximum size less than or equal to
500g, and a second gear, `other`, that catches everything else. The
position of the knife-edge for both gears will occur at 0.05 x the
maximum size i.e. the selectivity parameters will be different for each
species and will depend on the maximum size.

For this we will need to change the `gear_params` data frame. If we want
to keep the original model, we should first make a copy before making
modifications.

``` r
params_multi_gear <- params
```

To start with we need to know what the maximum sizes of the species in
the model are so we can determine the knife-edge positions for each
species. These are stored in the `w_max` column of the `species_params`
data frame inside the `params` object. We can then use these maximum
sizes to set a vector of knife edge sizes that are 0.05 times the
maximum size:

``` r
gear_params(params_multi_gear)$knife_edge_size <- 
  species_params(params)$w_max * 0.05
```

Now we want to assign each species to either the *industrial* or *other*
gear.

``` r
gear <- rep("Industrial", nrow(species_params(params)))
gear[species_params(params)$w_max > 500] <- "Other"
gear_params(params_multi_gear)$gear <- gear
```

To check what has just happened let us look at the new gear parameter
data frame. Its print method shows the four columns that every gear
needs and only names the selectivity parameters at the bottom, so we
print the knife-edge sizes separately:

``` r
gear_params(params_multi_gear)
```

    ## An object of class "gear_params" containing 10 gear-species pairs for 2 gears:
    ##        gear species   sel_func catchability
    ##  Industrial       1 knife_edge            1
    ##  Industrial       2 knife_edge            1
    ##  Industrial       3 knife_edge            1
    ##  Industrial       4 knife_edge            1
    ##       Other       5 knife_edge            1
    ##       Other       6 knife_edge            1
    ##       Other       7 knife_edge            1
    ##       Other       8 knife_edge            1
    ##       Other       9 knife_edge            1
    ##       Other      10 knife_edge            1
    ## With 1 other parameters: knife_edge_size

``` r
gear_params(params_multi_gear)$knife_edge_size
```

    ## 1, Industrial 2, Industrial 3, Industrial 4, Industrial      5, Other 
    ##     0.4456255     1.2559432     3.5397289     9.9763116    28.1170663 
    ##      6, Other      7, Other      8, Other      9, Other     10, Other 
    ##    79.2446596   223.3417961   629.4627059  1774.0669462  5000.0000000

Having created our `MizerParams` object with multiple gears, we can now
turn our attention to running a projection with multiple gears. In our
previous examples of calling
[`project()`](https://sizespectrum.org/mizer/reference/project.md) we
have specified the fishing effort with the `effort` argument using a
single value. This fixes the fishing effort for all gears in the model,
for all time steps. We can do this with our multi-gear parameter object:

``` r
sim_multi_gear <- project(params_multi_gear, t_max = 75, effort = 0.5)
```

By plotting this you can see that the fishing mortality for each species
now has a different selectivity pattern, and that the position of the
selectivity knife-edge is given by the maximum size of the species.

``` r
plot(sim_multi_gear)
```

![Summary plot of the trait-based model with multiple gears and a single
effort.](trait_model_files/figure-html/plot_multi_gear_trait_single_effort-1.png)

For the industrial fishery we said that we only wanted species with a
maximum size of 500 g or less to be fished. There are several ways of
specifying the `effort` argument for
[`project()`](https://sizespectrum.org/mizer/reference/project.md) .
Above we specified a single value that was used for all gears, for all
time steps. It is also possible to specify a separate effort for each
gear that will be used for all time steps. To do this we pass in effort
as a named vector. Here we set the effort for the *Industrial* gear to
0.75, and the effort of the *Other* gear to 0 (effectively turning it
off).

``` r
sim_multi_gear <- project(params_multi_gear, t_max = 75,
    effort = c(Industrial = 0.75, Other = 0))
```

Now you can see that the *Industrial* gear has been operating and that
fishing mortality for species larger than 500 g is 0.

``` r
plot(sim_multi_gear)
```

![Summary plot of the trait-based model with only the Industrial gear
active.](trait_model_files/figure-html/plot_multi_gear_trait-1.png)

## The impact of industrial fishing

In the previous section we set up and ran a model in which an industrial
fishery was operating that only selected smaller species. We can now
answer the question: what is the impact of such a fishery? We can again
compare abundances of the fished (`sim_industrial1`) and unfished
(`sim_industrial0`) cases:

``` r
sim_industrial0 <- project(params_multi_gear, t_max = 75, effort = 0)
sim_industrial1 <- project(params_multi_gear, t_max = 75,
    effort = c(Industrial = 0.75, Other = 0))
```

And compare them in the same way as before:

``` r
plotSpectraRelative(sim_industrial0, sim_industrial1, total = TRUE,
                    resource = FALSE, highlight = "Total")
```

![The relative difference between the industrially fished and the
unfished trait-based model at each size, showing a cascade that spreads
both upwards and downwards from the small species that are being caught.
Three of the four targeted species lie flat at -2, having gone
extinct.](trait_model_files/figure-html/plot_relative_comm_abund_industrial-1.png)

This shows another trophic cascade, although this time one driven by
fishing the small species at the lower end of the spectrum, not the
largest individuals as before. This trophic cascade acts in both
directions. The cascade upwards is driven by the lack of food for
predators leading to smaller realised maximum sizes. The cascade
downwards has the same mechanism as fishing on large fish, a combination
of predation mortality and food limitation.

Here too the fishing is hard enough to remove species outright: three of
the four species the *Industrial* gear catches go extinct, which is why
several of the coloured lines sit flat at \\-2\\. The larger species
that the gear does not touch are affected only through the food web:
their lines stay between roughly \\-1.6\\ and \\+1.2\\, well short of
the \\-2\\ that marks a disappearance.

The next section explains how to setup the more general [multispecies
model.](https://sizespectrum.org/mizer/articles/multispecies_model.md)
