# The Community Model

The simplest version of the size spectrum model is the community model.
In the community model, individuals are only characterized by their size
and are represented by a single group representing an across-species
aggregate. Reproduction is not considered; the reproduction rate \\R\\
is set to be constant. The resource spectrum only extends to the start
of the community spectrum. Metabolism is turned off.

In this section we describe how a community model can be set up and
projected through time. We then use a community model to illustrate the
idea of a trophic cascade. Due to the relative simplicity of this type
of model, it is useful for gently introducing some of the concepts
behind the mizer package. Consequently, this section should hopefully
serve as an introduction to using mizer.

## Setting up a community model

The first stage in implementing a model using `mizer` is to create an
object of class `MizerParams`. This class contains the model parameters
including the life-history parameters of the species in the model, the
fishing selectivity functions and the parameters of the resource
spectrum.

To avoid having to make a `MizerParams` object directly, the
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
wrapper function, has been provided that conveniently creates a
`MizerParams` object specifically for a community model. The
documentation for the function can be seen by clicking on the function
name anywhere it appears on this page:
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md).

As can be seen in the help page, the function can take many arguments.
We can ignore most of these for the moment as they almost all come with
default values.

The arguments that you should pay attention to are: `z0` (the external
mortality rate), `alpha` (the assimilation efficiency of the community),
`f0` (the average feeding level of the community, which is used to
calculate \\\gamma\\) and `h` (the coefficient of the maximum intake
rate).

Although default values for these parameters are provided, you are
encouraged to explore how changing the values affects the simulated
community. For example, the default value of `z0` is \\0.1\\. Increasing
this value effectively ‘shortens’ the length of the community spectrum.

The
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
function is called by passing in the arguments by name. Any parameter
that is not passed in is set to the default value. For example, the
following line sets up the parameters with `z0` = 0.05, `f0` = 0.5. All
other parameters will have their default value:

``` r

params <- newCommunityParams(z0 = 0.05, f0 = 0.5)
```

Calling the function creates and returns an object of type
`MizerParams`. We can check this using the
[`class()`](https://rdrr.io/r/base/class.html) function.

``` r

class(params)
```

    ## [1] "MizerParams"

To work with mizer you do not have to worry about how this object is
realised internally. All you need to know is that it contains all the
information needed to run model simulations and that mizer provides you
with functions to access, use and modify this information.

A very brief description of the model contained in the MizerParams
object can be seen by calling the
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md)
method on it:

``` r

summary(params)
```

    ## An object of class "MizerParams" 
    ## mizer version: 3.4.0
    ## Created: 2026-09-02 10:02:30
    ## Modified: 2026-09-02 10:02:30
    ## Consumer size spectrum:
    ##  minimum size:   0.001
    ##  maximum size:   1e+06
    ##  no. size bins:  100
    ## Resource size spectrum:
    ##  minimum size:   8.11131e-11
    ##  maximum size:   0.000811131
    ##  no. size bins:  78  (178 size bins in total)
    ## Steady state:
    ##  biomass drift:  4.2 /year   (not at steady state, largest in the resource - run tuneSteadyState())
    ## Species details:
    ## An object of class "species_params" containing parameters for 1 species:
    ##    species w_inf  w_mat w_min  f0 beta sigma
    ##  Community 1e+06 250000 0.001 0.5  100     2
    ## 
    ## Fishing gear details:
    ## Gear          Effort  Target species 
    ##  ----------------------------------
    ## Community      0.00   Community

In the summary you can see that the size range of the community spectrum
has been set from \\0.001\\ to \\10^{6}\\ and this is divided into
\\100\\ size bins. Similar information is available for the resource
spectrum. Additionally, the community is made up of only one species,
called *Community*, which has an maximum size of \\10^{6}\\ and a
preferred predator prey mass ratio of \\100\\. The `w_mat` parameter is
filled in with its usual default but plays no role in a community model:
the reproduction rate is held constant and individuals invest none of
their income into reproduction. These values have all been set by
default using the
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
function.

The summary also reports how far the model is from steady state,
measured by the rate at which the biomass is currently drifting. We have
not tuned this model to a steady state, so the abundances will change
when we start to project the model through time. For a model you intend
to use for projections you would normally first find its steady state;
see the guide on [calibrating a
model](https://sizespectrum.org/mizer/articles/guide-calibrate-model.md).

## Running the community model

The
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
function has given us a `MizerParams` object that contains all the
information mizer needs about the model community. We can use this to
perform a simulation and project the community through time. In the
`mizer` package, projections are performed using the
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
function. You can see the help page for
[`project()`](https://sizespectrum.org/mizer/reference/project.md) for
more details and it is described fully in [the guide on running a
simulation.](https://sizespectrum.org/mizer/articles/guide-run-simulation.md)
We will ignore the details for the moment and just use
[`project()`](https://sizespectrum.org/mizer/reference/project.md) to
run some simple projections. The arguments for
[`project()`](https://sizespectrum.org/mizer/reference/project.md) that
we need to be concerned with are `effort`, which determines the fishing
effort (and therefore fishing mortality) through time, and `t_max`,
which is the number of years we want to project into the future. Initial
population abundances have been set automatically (by the
[`get_initial_n()`](https://sizespectrum.org/mizer/reference/get_initial_n.md)
function). It is possible to set your own initial abundances but we will
not do this here.

To run a projection for 50 years, with no fishing effort (i.e. we want
to model an unexploited community) we run:

``` r

sim <- project(params, t_max = 50, effort = 0)
```

The resulting object, `sim`, is of type `MizerSim`.

``` r

class(sim)
```

    ## [1] "MizerSim"

This class holds the results of the simulation, including the community
and resource abundances at size through time, as well as the original
model parameters. It is explained in detail in [the guide on running a
simulation.](https://sizespectrum.org/mizer/articles/guide-run-simulation.md)

After running the projection, it is possible to explore the results
using a range of plots and analyses. These are described fully in [the
guide on analysing and plotting
results.](https://sizespectrum.org/mizer/articles/guide-analyse-and-plot.md)

To quickly look at the results of the projection you can call the
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) method.
This plots the feeding level, predation mortality, fishing mortality and
abundance by size in the last time step of the simulation, and the total
biomass through time. Each of the plots can be shown individually if
desired.

``` r

plot(sim)
```

![The plot shows the results of the community model projection. The top
panel shows the total biomass of the community through time. The other
panels show the feeding level, predation mortality, fishing mortality
and abundance by size in the last time step of the
simulation.](community_model_files/figure-html/print_plot_comm_sim-1.png)

In the above plot there are several things going on that are worth
talking about. Looking at the total biomass of the community against
time, you can see that the biomass quickly reaches a stable equilibrium.
The other panels show what is happening at the last time step in the
simulation, which in this case is when the community is at equilibrium.
Fishing mortality is 0 because we set the `effort` argument to 0 when
running the simulation. The predation mortality rate is clearly a
function of size, with the smallest sizes experiencing the highest
levels of predation. The feeding level describes how satiated an
individual is, with 0 being unfed, and 1 being full satiated. The
feeding level at size will be strongly affected by the values of the
`f0` and `alpha` arguments passed to the
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
function.

The resource and community spectra are shown in the bottom panel of the
plot (the plotted resource spectrum has been truncated to make for a
better plot, but really extends all the way back to \\8.11\times
10^{-11}\\ g). You can see that the community spectrum forms a continuum
with the resource spectrum. This is strongly affected by the level of
fixed reproduction rate (the `reproduction` argument passed to
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md))

Note the hump in the biomass at the largest end of the community
spectrum. This is because the size spectrum model can be broadly
described as ‘big things eating little things’. Given this, what is
eating the very biggest things? Without fishing pressure, the mortality
of the largest individuals is only from the external mortality
(determined by the `z0` argument) and the mortality from predation is
almost 0. This is difficult to see in the summary plot because the
predation mortality is so high for the smaller individuals.

Each panel of the summary plot is also available as a plot function of
its own. Here we use
[`plotPredMort()`](https://sizespectrum.org/mizer/reference/plotPredMort.md)
to look at the predation mortality alone. Because the predation
mortality spans many orders of magnitude we ask for a logarithmic y axis
as well as the logarithmic x axis that these plots use by default:

``` r

plotPredMort(sim, log_y = TRUE)
```

![The plot shows the predation mortality at size in the final time step
of the community model projection. Both axes are on a log scale. The
predation mortality declines to almost zero for the largest
sizes.](community_model_files/figure-html/print_plot_comm_m2-1.png)

Now the decline of the predation mortality to almost zero for the
largest sizes is clearly visible. There are analogous functions for the
other panels of the summary plot:
[`plotFeedingLevel()`](https://sizespectrum.org/mizer/reference/plotFeedingLevel.md),
[`plotFMort()`](https://sizespectrum.org/mizer/reference/plotFMort.md),
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)
and
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md).
They all accept a `MizerSim` object, in which case they show the final
time step of the simulation, or a `MizerParams` object, in which case
they show the state that the `MizerParams` object describes.

If you want the numbers behind such a plot rather than the picture,
there is a `get...()` function for each of them. Here we use
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md):

``` r

pred_mort <- getPredMort(sim)
dim(pred_mort)
```

    ## [1]  51 100

The resulting `pred_mort` object contains the predation mortality at
time by species by size. Here we have only one species, so the species
dimension is dropped, leaving us with a two-dimensional array of time by
size. We projected the model for \\50\\ time steps but the length of the
time dimension is \\51\\ because the initial state is also stored.

Rather than picking the final row out of that array by hand, you extract
the state of the model at the final time step with
[`finalParams()`](https://sizespectrum.org/mizer/reference/getParams.md).
This returns a `MizerParams` object holding the abundances at the end of
the simulation, which you can then hand to any of the `get...()` or
`plot...()` functions:

``` r

pred_mort_final <- getPredMort(finalParams(sim))
```

The related function `getParams(sim, time_range = ...)` gives you the
state averaged over any period of the simulation, and
`initialParams(sim)` gives the state the simulation started from.

The object returned by
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md)
is a mizer array that knows its own units and the model it came from, so
it can plot itself. The following produces the same plot as
[`plotPredMort()`](https://sizespectrum.org/mizer/reference/plotPredMort.md)
above:

``` r

plot(pred_mort_final, log_y = TRUE)
```

![The predation mortality at size in the final time step, plotted
directly from the array returned by
getPredMort().](community_model_files/figure-html/plot_m2_array-1.png)

This works for the result of any of the summary functions, so you rarely
need to write plotting code of your own. See the guide on [analysing and
plotting
results](https://sizespectrum.org/mizer/articles/guide-analyse-and-plot.md)
for the full set of functions and the arguments they share.

## Example of a trophic cascade with the community model

It is possible to use the community model to simulate a trophic cascade.
To do this we need to perform two simulations, one with fishing and one
without.

This means we need to consider how fishing is handled in `mizer`. The
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
function automatically sets the fishing selectivity to have a knife-edge
shape, with only individuals larger than 1 kg selected (the size at the
knife-edge can be changed by setting the `knife_edge_size` argument).
Although it is possible to change the selectivity function, here we will
use the default knife-edge selectivity. We set up the parameter object
with default parameters:

``` r

params_knife <- newCommunityParams()
```

First we perform a simulation without fishing in the same way we did
above by setting the `effort` argument to 0:

``` r

sim0 <- project(params_knife, effort = 0, t_max = 50)
```

Now we want to simulate again, this time with fishing. In the
simulations, fishing mortality is calculated as the product of the
fishing selectivity, effort and catchability (see [the guide on setting
up
fishing](https://sizespectrum.org/mizer/articles/guide-set-up-fishing.md)
for more details). By default catchability is set to 1. This means that
a fishing effort of 1 will result in a fishing mortality of 1/year for
fully selected sizes. Here we run a simulation with fishing effort set
to 1 for the duration of the simulation:

``` r

sim1 <- project(params_knife, effort = 1, t_max = 50)
```

You can compare the difference between these scenarios by using the
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) method as
before. Of particular interest is the fishing mortality at size. The
knife-edge selectivity at 1000 g can be clearly seen and an effort of 1
has resulted in a fishing mortality of 1 for the fully selected sizes.

``` r

plot(sim1, biomass = TRUE, per_log_size = TRUE)
```

![A summary plot of the simulation with
fishing.](community_model_files/figure-html/print_plot_comm_fmort-1.png)

To explore the presence of a trophic cascade, we are interested in the
change in abundance when the community is fished compared to when it is
not fished. The abundances at size are held in the `MizerSim` object and
are obtained with the
[`N()`](https://sizespectrum.org/mizer/reference/N.md) function, which
returns a three dimensional array with dimensions time x species x size.
Here we have 51 time steps (50 from the simulation plus one which stores
the initial state), 1 species and 100 sizes:

``` r

dim(N(sim0))
```

    ## [1]  51   1 100

We are interested in the abundances at the final time step. Instead of
working out the index of the final time yourself, use
[`finalN()`](https://sizespectrum.org/mizer/reference/finalN.md), which
is the counterpart of
[`finalParams()`](https://sizespectrum.org/mizer/reference/getParams.md)
for the abundances alone. The ratio of the two gives us the abundance in
the fished community relative to the unfished one:

``` r

relative_abundance <- finalN(sim1) / finalN(sim0)
range(relative_abundance)
```

    ## [1] 3.496447e-21 2.548582e+00

We do not need to plot that by hand, though, because mizer has a
function for exactly this comparison.
[`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md)
takes two models or simulations and plots the difference between their
spectra relative to their average, i.e. \\2(N_1(w) - N_0(w)) / (N_1(w) +
N_0(w))\\. This is a more comfortable measure than the plain ratio,
because it stays in the range from \\-2\\ to \\2\\ however extreme the
change is: values above zero mark sizes that are more abundant in the
fished community, values below zero sizes that are less abundant.
Because the factors of \\w\\ cancel in a relative difference, it makes
no difference here whether you think in terms of number density or
biomass density.

``` r

plotSpectraRelative(sim0, sim1)
```

![The plot shows the relative difference between the abundance of the
fished and the unfished community at each size, with a strong reduction
above 1000 g and alternating increases and decreases at smaller
sizes.](community_model_files/figure-html/plot_relative_comm_abund-1.png)

If you would rather see the two spectra themselves instead of their
difference,
[`plotSpectra2()`](https://sizespectrum.org/mizer/reference/plotSpectra2.md)
overlays them. Fishing removes the large individuals almost entirely, so
we use `ylim` to cut off the very bottom of the y axis and keep the
interesting part of the plot readable:

``` r

plotSpectra2(sim0, sim1, name1 = "Unfished", name2 = "Fished",
             per_log_size = TRUE, ylim = c(1, NA))
```

![The biomass density spectra of the unfished and the fished community
overlaid, showing the collapse of the fished spectrum above 1000 g and
the alternating bumps and dips it develops at smaller
sizes.](community_model_files/figure-html/plot_spectra2_comm-1.png)

The impact of fishing on individuals larger than 1000 g can be clearly
seen. The fishing pressure lowers the abundance of large fish. This then
relieves the predation pressure on their smaller prey (the preferred
predator-prey size ratio is given by the \\\beta\\ parameter, which is
set to 100 by default), leading to an increase in their abundance. This
in turn increases the predation mortality on *their* smaller prey, which
reduces their abundance and so on down the spectrum.

## The impact of changing \\\sigma\\

As described above, the \\\sigma\\ parameter determines the width of the
predator prey size preference. Here we take a look at how changing the
value of \\\sigma\\ can affect the dynamics of the community.

In the examples above, \\\sigma\\ is set in the
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
function by default to a value of \\2\\. We can see this by looking at
the `sigma` column of the species parameter data frame that is contained
in the MizerParams object:

``` r

species_params(params)$sigma
```

    ## Community 
    ##         2

When projected through time, the community abundances converge to a
stable equilibrium. What happens if we reduce the value of \\\sigma\\,
for example by setting it to 1.0? We can do this by passing in the new
value of \\\sigma\\ into
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md).

``` r

params_sigma1 <- newCommunityParams(sigma = 1)
```

We want to project this new model through time using the
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
function. Here we project the new parameter object for 50 time steps
without fishing and save at intervals of 0.1 years (`t_save = 0.1`):

``` r

sim_sigma1 <- project(params_sigma1, effort = 0, t_max = 50, 
                      dt = 0.01, t_save = 0.1)
```

Note that we have introduced a new argument: \\dt\\. This is the step
size of the solver. It does not have anything to do with the biology in
the model. It only affects the internal engine of
[`project()`](https://sizespectrum.org/mizer/reference/project.md) that
performs the projection. As you can see in the underlying model
equations in [the model description
section](https://sizespectrum.org/mizer/articles/model_description.md),
the model is formulated in continuous time. Therefore, to project it
forward,
[`project()`](https://sizespectrum.org/mizer/reference/project.md) must
solve the system of equations using numerical methods. The quality of
these methods is strongly affected by \\dt\\. The default value of
\\dt\\ is 0.1, which will be fine for most of the projections we run in
this document. Here it is necessary to reduce the value to 0.01 to avoid
introducing any artefacts into the projected values. Decreasing \\dt\\
increases the time it takes to run a projection.

Let’s take a look at how the abundances change through time. We can do
this with the
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md)
function:

``` r

plotBiomass(sim_sigma1)
```

![The plot shows the biomass of the community through time when the
value of sigma is set to 1.0. The abundances no longer converge to a
stable equilibrium and the dynamics appear to be
chaotic.](community_model_files/figure-html/plot_comm_biomass_sigma1-1.png)

The plot above shows that abundances of the community no longer converge
to a stable equilibrium and the dynamics appear to be chaotic. The
ecological significance of the change in dynamics, and of the ability of
simple community models to show chaotic behaviour, is still being
debated. It can be argued that the size of the oscillations are too
large to be ‘true’. Additionally, when a trait-based model is
implemented, the magnitude of the oscillations are much smaller. Mizer
provides tools for investigating when a model loses its stability in
this way; see the guide on [analysing the stability of a
model](https://sizespectrum.org/mizer/articles/guide-analyse-stability.md).

The next section is about [the trait based
model.](https://sizespectrum.org/mizer/articles/trait_model.md)
