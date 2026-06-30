# The Multi Species Model

The previous sections have used wrapper functions to set up MizerParams
objects that are appropriate for single-species, community- and
trait-based models. We now turn our attention to multispecies, or
species-specific, models. These are potentially more complicated than
the community and trait-based models and use the full power of the
`mizer` package.

In multispecies type models multiple species are resolved. However,
unlike in the trait-based model which also resolves multiple species,
the species parameters will be those of real-world species. There are
several advantages to this approach. As well as investigating the
community as a whole (as was done for the community and trait-based
models), we are able to investigate the dynamics of individual species.
This means that species specific management rules can be tested and
species specific metrics, such as yield, can be compared to reference
levels.

A multispecies model can take more effort to set up. For example, each
species will have different life-history parameters; there may be
multiple gear types with different selectivities targeting different
groups of species; the fishing effort of each gear may change with time
instead of just being constant (which has been the case in the
simulations we have looked at so far); the interactions between the
species needs to be considered.

In later sections we build up a multispecies model for the North Sea. To
effectively use `mizer` for a multispecies model we are going to have to
take a closer look at the `MizerParams` class and the
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
function. This will all be done in the context of examples so hopefully
everything will be clear.

We also take a closer look at some of the summary plots and analyses
that can be performed, for example, calculating a range of size-based
indicators.

## Setting up a multispecies model

### Overview

The MizerParams class is used for storing model parameters. We have
already met the MizerParams class when we looked at community and
trait-based models. However, to set up a multispecies model we will need
to specify many more parameters.This is probably the most complicated
part of using the mizer package, so we will take it slowly.

A `MizerParams` object stores the:

- life-history parameters of the species in the community, such as
  maximum size \\w\_{max}\\;
- size-based biological parameters for the species, such as the search
  volume;
- density-dependent reproduction functions and parameters of each
  species;
- an interaction matrix to describe the spatial overlap of pairs of
  species;
- parameters relating to the growth and dynamics of the resource
  spectrum;
- fishing gear parameters: selectivity and catchability.

Note that the `MizerParams` class does not store any parameters that can
vary through time, such as fishing effort or population abundance. These
are stored in the `MizerSim` class which we will come to later in [the
section on running a
simulation.](https://sizespectrum.org/mizer/articles/running_a_simulation.html#sec:projection)

Although the `MizerParams` class contains a lot of information, it is
relatively straightforward to set up and use. Objects of class
`MizerParams` are created using the constructor method
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
(this constructor method was called MizerParams() in previous version of
mizer). This constructor method can take many arguments. However,
creation is simplified because many of the arguments have default
values.

In the rest of this section we look at the main arguments to the
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
function. To help understand how the constructor is used and how the
`MizerParams` class relates to the equations given in [the model
description
section,](https://sizespectrum.org/mizer/articles/model_description.md)
there is an example section where we create example parameter objects
using data that comes with the `mizer` package.

### The species parameters

Although many of the arguments used when creating a `MizerParams` object
are optional, there is one argument that must be supplied by the user:
the *species specific parameters*. These are stored in a single
`data.frame` object. The `data.frame` is arranged species by parameter,
so each column is a parameter and each row has the parameters for one of
the species in the model. Although it is possible to create the
data.frame by hand in R, it is probably easier to create the data
externally as a .csv file (perhaps using a suitable open source
spreadsheet such as LibreOffice) and then read the data into R.

For each species in the model community there are certain parameters
that are essential and that do not have default values. The user must
provide values for these parameters. There are also some essential
parameters that have default values, such as the selectivity function
parameters, and some that are calculated internally using default
relationships if not explicitly provided. These defaults are used if the
parameters are not found in the data.frame.

The essential columns of the species parameters data.frame that have no
default values are: `species`, the names of the species in the community
and `w_inf`, the von Bertalanffy asymptotic mass of the species. (The
computational upper size boundary `w_max` is not essential and defaults
to `1.5 * w_inf`.)

### The gear parameters

In `mizer`, fishing mortality is imposed on species by fishing gears.
The total fishing mortality is obtained by summing over the mortality
from all gears, \\\begin{equation} % {#eq:muf} \mu\_{f.i}(w) = \sum_g
F\_{g,i}(w), \end{equation}\\ where the fishing mortality
\\F\_{g,i}(w)\\ imposed by gear \\g\\ on species \\i\\ at size \\w\\ is
calculated as: \\\begin{equation} % {#eq:sel} F\_{g,i}(w) = S\_{g,i}(w)
Q\_{g,i} E\_{g} \end{equation}\\ where \\S\\ is the selectivity by
species, gear and size, \\Q\\ is the catchability by species and gear
and \\E\\ is the fishing effort by gear. The selectivity at size has a
range between 0 (not selected at that size) to 1 (fully selected at that
size). Catchability is used as an additional scalar to make the link
between gear selectivity, fishing effort and fishing mortality. For
example, it can be set so that an effort of 1 gives a desired fishing
mortality. In this way effort can then be specified relative to a ‘base
effort’, e.g. the effort in a particular year.

Selectivity and catchability are stored as arrays in the MizerParams
object. However the user does not have to create these arrays by hand if
they provide a data frame with the necessary information. In particular
the selectivity can be calculate by specifying functions for the
selectivity curves. Mizer provides a range of such selectivity functions
and the user just needs to specify their parameters for each gear and
each species in the `gear_params` data frame. All the details can be
found on the help page for
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md).

Fishing effort is not stored in the MizerParams object. Instead, effort
is set when the simulation is run and can vary through time (see [the
section on running a
simulation](https://sizespectrum.org/mizer/articles/running_a_simulation.md)).

### Example of making `MizerParams` objects

As mentioned in the preceding sections, an object of `MizerParams` is
created by using the
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
constructor method.

The first step is to prepare the species specific parameter data.frame.
As mentioned above, one way of doing this is to use a spreadsheet and
save it as a .csv file. We will use this approach here. An example .csv
file has been included in the package. This contains the species
parameters for a multispecies North Sea model. The location of the file
can be found by running

``` r

params_location <- system.file("extdata", "NS_species_params.csv",
                               package = "mizer")
```

This file can be opened with most spreadsheets or a text editor for you
to inspect. This can be loaded into R with

``` r

species_params <- read.csv(params_location)
```

This reads the .csv file into R in the form of a data.frame. You can
check this with the `class`:

``` r

class(species_params)
```

    ## [1] "data.frame"

Let’s have a look at the data frame:

``` r

species_params
```

    ##    species w_inf w_mat   beta sigma    R_max  k_vb
    ## 1    Sprat    33    13  51076   0.8 7.38e+11 0.681
    ## 2  Sandeel    36     4 398849   1.9 4.10e+11 1.000
    ## 3   N.pout   100    23     22   1.5 1.05e+13 0.849
    ## 4  Herring   334    99 280540   3.2 1.11e+12 0.606
    ## 5      Dab   324    21    191   1.9 1.12e+10 0.536
    ## 6  Whiting  1192    75     22   1.5 5.48e+11 0.323
    ## 7     Sole   866    78    381   1.9 3.87e+10 0.284
    ## 8  Gurnard   668    39    283   1.8 1.65e+12 0.266
    ## 9   Plaice  2976   105    113   1.6 4.08e+14 0.122
    ## 10 Haddock  3485   165    558   2.1 1.84e+12 0.271
    ## 11     Cod 40044  1606     66   1.3 8.26e+09 0.216
    ## 12  Saithe 16856  1076     40   1.1 1.12e+11 0.175

You can see that there are \\12\\ species and \\7\\ columns of
parameters: `species`, `w_inf`,`w_mat`,`beta`,`sigma`,`R_max` and
`k_vb`.

Of these parameters, `species` and `w_inf` are essential and have no
default values (as described in [the section on species
parameters](https://sizespectrum.org/mizer/articles/multispecies_model.html#sec:species_parameters_dataframe)).
`w_inf` is the von Bertalanffy asymptotic size of the species, `w_mat`
is its maturity size, and `beta` and `sigma` are parameters of the
predation kernel (by default mizer uses a log-normal predation kernel).
`R_max` is a parameter introducing additional density dependence into
reproduction parameter using a Beverton-Holt type function (see
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
for details). The final column, `k_vb`, will be used to calculate values
for `h` and then `gamma`. This column is only essential here because the
`h` and `gamma` are not included in the data.frame. It would also have
been possible to include `h` and `gamma` columns in the data.frame and
not include the `k_vb` column.

The values of the non-essential species specific parameters, like for
example `alpha`, `k`, `ks`, `z0`, `w_min` and `erepro`, were not
included in the data.frame. This means that the default values will be
automatically used when we create the `MizerParams` object.

For this example we will not set up gear parameters. There are no
columns describing the fishing selectivity. There is no `sel_func`
column to determine the selectivity function. This means that the
default selectivity function, `knife_edge`, will be used. As mentioned
in [the section on fishing
gears](https://sizespectrum.org/mizer/articles/multispecies_model.html#sec:fishing_gear),
this function also needs another argument, `knife_edge_size`. This is
not present in the data.frame and so it will be set to the default value
of `w_mat`. Also, there is no `catchability` column so a default value
for `catchability` of 1 will be used for all gears and species.

To create the `MizerParams` object we pass the species parameter
data.frame into the
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
constructor method:

``` r

params <- newMultispeciesParams(species_params)
```

    ## Because you have n != p, the default value for `h` is not very good.
    ## Because the age at maturity is not known, I need to fall back to using
    ## von Bertalanffy parameters, where available, and this is not reliable.
    ## No ks column so calculating from critical feeding level.
    ## Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.
    ## Using f0, h, lambda, kappa and the predation kernel to calculate gamma.

We have just created a `MizerParams` object:

``` r

class(params)
```

    ## [1] "MizerParams"
    ## attr(,"package")
    ## [1] "mizer"

The MizerParams object also stores a copy of the species parameter data
frame that we provided. We can look at it with
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md):

``` r

species_params(params)
```

    ##         species w_inf w_mat   beta sigma    R_max  k_vb         n   p w_min
    ## Sprat     Sprat    33    13  51076   0.8 7.38e+11 0.681 0.6666667 0.7 0.001
    ## Sandeel Sandeel    36     4 398849   1.9 4.10e+11 1.000 0.6666667 0.7 0.001
    ## N.pout   N.pout   100    23     22   1.5 1.05e+13 0.849 0.6666667 0.7 0.001
    ## Herring Herring   334    99 280540   3.2 1.11e+12 0.606 0.6666667 0.7 0.001
    ## Dab         Dab   324    21    191   1.9 1.12e+10 0.536 0.6666667 0.7 0.001
    ## Whiting Whiting  1192    75     22   1.5 5.48e+11 0.323 0.6666667 0.7 0.001
    ## Sole       Sole   866    78    381   1.9 3.87e+10 0.284 0.6666667 0.7 0.001
    ## Gurnard Gurnard   668    39    283   1.8 1.65e+12 0.266 0.6666667 0.7 0.001
    ## Plaice   Plaice  2976   105    113   1.6 4.08e+14 0.122 0.6666667 0.7 0.001
    ## Haddock Haddock  3485   165    558   2.1 1.84e+12 0.271 0.6666667 0.7 0.001
    ## Cod         Cod 40044  1606     66   1.3 8.26e+09 0.216 0.6666667 0.7 0.001
    ## Saithe   Saithe 16856  1076     40   1.1 1.12e+11 0.175 0.6666667 0.7 0.001
    ##           w_max w_repro_max alpha interaction_resource z_ext          d E_ext
    ## Sprat      49.5          33   0.6                    1     0 -0.3333333     0
    ## Sandeel    54.0          36   0.6                    1     0 -0.3333333     0
    ## N.pout    150.0         100   0.6                    1     0 -0.3333333     0
    ## Herring   501.0         334   0.6                    1     0 -0.3333333     0
    ## Dab       486.0         324   0.6                    1     0 -0.3333333     0
    ## Whiting  1788.0        1192   0.6                    1     0 -0.3333333     0
    ## Sole     1299.0         866   0.6                    1     0 -0.3333333     0
    ## Gurnard  1002.0         668   0.6                    1     0 -0.3333333     0
    ## Plaice   4464.0        2976   0.6                    1     0 -0.3333333     0
    ## Haddock  5227.5        3485   0.6                    1     0 -0.3333333     0
    ## Cod     60066.0       40044   0.6                    1     0 -0.3333333     0
    ## Saithe  25284.0       16856   0.6                    1     0 -0.3333333     0
    ##         D_ext is_background pred_kernel_type        h k       ks         z0
    ## Sprat       0         FALSE        lognormal 14.51026 0 1.598545 0.18705957
    ## Sandeel     0         FALSE        lognormal 28.36951 0 3.250607 0.18171206
    ## N.pout      0         FALSE        lognormal 30.69918 0 3.318311 0.12926608
    ## Herring     0         FALSE        lognormal 31.20041 0 3.212332 0.08647736
    ## Dab         0         FALSE        lognormal 34.68295 0 3.760307 0.08735805
    ## Whiting     0         FALSE        lognormal 32.78322 0 3.406676 0.05658819
    ## Sole        0         FALSE        lognormal 24.90951 0 2.585095 0.06294752
    ## Gurnard     0         FALSE        lognormal 22.29126 0 2.367448 0.06863713
    ## Plaice      0         FALSE        lognormal 17.71691 0 1.820523 0.04171321
    ## Haddock     0         FALSE        lognormal 40.62144 0 4.111691 0.03957464
    ## Cod         0         FALSE        lognormal 74.81794 0 7.019866 0.01753768
    ## Saithe      0         FALSE        lognormal 43.50194 0 4.136466 0.02340093
    ##                 q        gamma     w_mat25 m erepro
    ## Sprat   0.7166667 5.753738e-11   11.647460 1      1
    ## Sandeel 0.7166667 4.258153e-11    3.583834 1      1
    ## N.pout  0.7166667 9.729943e-11   20.607045 1      1
    ## Herring 0.7166667 2.806634e-11   88.699888 1      1
    ## Dab     0.7166667 7.647896e-11   18.815128 1      1
    ## Whiting 0.7166667 1.039047e-10   67.196884 1      1
    ## Sole    0.7166667 5.297276e-11   69.884760 1      1
    ## Gurnard 0.7166667 5.081126e-11   34.942380 1      1
    ## Plaice  0.7166667 4.764030e-11   94.075638 1      1
    ## Haddock 0.7166667 7.662872e-11  147.833146 1      1
    ## Cod     0.7166667 2.544302e-10 1438.909287 1      1
    ## Saithe  0.7166667 1.793361e-10  964.051303 1      1

We can see that this returns the original species data.frame (with
`w_inf` and so on), plus any default values that may not have been
included in the original data.frame. For example, we can see that there
are now columns for `alpha` and `h` and `gamma` etc.

Also note how the default fishing gears have been set up. Even though we
did not provide a gear parameter data frame, the MizerParams object has
one that we can access with

``` r

gear_params(params)
```

    ##                          species            gear   sel_func knife_edge_size
    ## Sprat, knife_edge_gear     Sprat knife_edge_gear knife_edge              13
    ## Sandeel, knife_edge_gear Sandeel knife_edge_gear knife_edge               4
    ## N.pout, knife_edge_gear   N.pout knife_edge_gear knife_edge              23
    ## Herring, knife_edge_gear Herring knife_edge_gear knife_edge              99
    ## Dab, knife_edge_gear         Dab knife_edge_gear knife_edge              21
    ## Whiting, knife_edge_gear Whiting knife_edge_gear knife_edge              75
    ## Sole, knife_edge_gear       Sole knife_edge_gear knife_edge              78
    ## Gurnard, knife_edge_gear Gurnard knife_edge_gear knife_edge              39
    ## Plaice, knife_edge_gear   Plaice knife_edge_gear knife_edge             105
    ## Haddock, knife_edge_gear Haddock knife_edge_gear knife_edge             165
    ## Cod, knife_edge_gear         Cod knife_edge_gear knife_edge            1606
    ## Saithe, knife_edge_gear   Saithe knife_edge_gear knife_edge            1076
    ##                          catchability
    ## Sprat, knife_edge_gear              1
    ## Sandeel, knife_edge_gear            1
    ## N.pout, knife_edge_gear             1
    ## Herring, knife_edge_gear            1
    ## Dab, knife_edge_gear                1
    ## Whiting, knife_edge_gear            1
    ## Sole, knife_edge_gear               1
    ## Gurnard, knife_edge_gear            1
    ## Plaice, knife_edge_gear             1
    ## Haddock, knife_edge_gear            1
    ## Cod, knife_edge_gear                1
    ## Saithe, knife_edge_gear             1

All species are caught by a gear called “knife_edge_gear”. The
selectivity function for each fishing gear has been set in the
`sel_func` column to the default function,
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md).
A `catchability` column has been added with a default value of 1 for
each of the species that the gear catches. An example of setting the
catchability by hand can be seen in [the section on the North
Sea.](https://sizespectrum.org/mizer/articles/a_multispecies_model_of_the_north_sea.md)

There is a
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md)
method for `MizerParams` objects which prints a useful summary of the
model parameters:

``` r

summary(params)
```

    ## An object of class "MizerParams" 
    ## mizer version: 3.1.0
    ## Created: 2026-06-25 21:02:52
    ## Modified: 2026-06-25 21:02:52
    ## Consumer size spectrum:
    ##  minimum size:   0.001
    ##  maximum size:   60066
    ##  no. size bins:  100
    ## Resource size spectrum:
    ##  minimum size:   9.20914e-13
    ##  maximum size:   8.48399
    ##  no. size bins:  166 (215 size bins in total)
    ## Species details:
    ##    species w_inf w_mat w_min   beta sigma
    ## 1    Sprat    33    13 0.001  51076   0.8
    ## 2  Sandeel    36     4 0.001 398849   1.9
    ## 3   N.pout   100    23 0.001     22   1.5
    ## 4  Herring   334    99 0.001 280540   3.2
    ## 5      Dab   324    21 0.001    191   1.9
    ## 6  Whiting  1192    75 0.001     22   1.5
    ## 7     Sole   866    78 0.001    381   1.9
    ## 8  Gurnard   668    39 0.001    283   1.8
    ## 9   Plaice  2976   105 0.001    113   1.6
    ## 10 Haddock  3485   165 0.001    558   2.1
    ## 11     Cod 40044  1606 0.001     66   1.3
    ## 12  Saithe 16856  1076 0.001     40   1.1
    ## 
    ## Fishing gear details:
    ## Gear          Effort  Target species 
    ##  ----------------------------------
    ## knife_edge_gear 0.00   Sprat, Sandeel, N.pout, Herring, Dab, Whiting, Sole, Gurnard, Plaice, Haddock, Cod, Saithe

As well as giving a summary of the species in the model and what gear is
fishing what species, it gives a summary of the size structure of the
community. For example there are \\100\\ size classes in the community,
ranging from \\0.001\\ to \\6.01\times 10^{4}\\ . These values are
controlled by the arguments `no_w`, `min_w` and `max_w` respectively.
For example, if we wanted 200 size classes in the model we would use:

``` r

params200 <- newMultispeciesParams(species_params, no_w = 200)
summary(params200)
```

### Setting the interaction matrix

So far we have created a `MizerParams` object by passing in only the
species parameter data.frame argument. We did not specify an interaction
matrix. The interaction matrix describes the interaction of each pair of
species in the model. This can be viewed as a proxy for spatial
interaction e.g. to model predator-prey interaction that is not size
based. The values in the interaction matrix are used to scale the
encountered food in \[getEncounter()\] and the predation mortality rate
in \[getPredMort()\] (see [the section on predator-prey encounter
rate](https://sizespectrum.org/mizer/articles/model_description.html#sec:pref)
and on [predation
mortality](https://sizespectrum.org/mizer/articles/model_description.html#mortality)).

The entries of the interaction matrix are dimensionless numbers taking
values are between 0 (species do not overlap and therefore do not
interact with each other) to 1 (species overlap perfectly). By default
mizer sets all values to 1, implying that all species fully interact
with each other, i.e. the species are spread homogeneously across the
model area.

``` r

getInteraction(params)
```

    ## Warning: `getInteraction()` was deprecated in mizer 2.4.0.
    ## ℹ Please use `interaction_matrix()` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ##          prey
    ## predator  Sprat Sandeel N.pout Herring Dab Whiting Sole Gurnard Plaice Haddock
    ##   Sprat       1       1      1       1   1       1    1       1      1       1
    ##   Sandeel     1       1      1       1   1       1    1       1      1       1
    ##   N.pout      1       1      1       1   1       1    1       1      1       1
    ##   Herring     1       1      1       1   1       1    1       1      1       1
    ##   Dab         1       1      1       1   1       1    1       1      1       1
    ##   Whiting     1       1      1       1   1       1    1       1      1       1
    ##   Sole        1       1      1       1   1       1    1       1      1       1
    ##   Gurnard     1       1      1       1   1       1    1       1      1       1
    ##   Plaice      1       1      1       1   1       1    1       1      1       1
    ##   Haddock     1       1      1       1   1       1    1       1      1       1
    ##   Cod         1       1      1       1   1       1    1       1      1       1
    ##   Saithe      1       1      1       1   1       1    1       1      1       1
    ##          prey
    ## predator  Cod Saithe
    ##   Sprat     1      1
    ##   Sandeel   1      1
    ##   N.pout    1      1
    ##   Herring   1      1
    ##   Dab       1      1
    ##   Whiting   1      1
    ##   Sole      1      1
    ##   Gurnard   1      1
    ##   Plaice    1      1
    ##   Haddock   1      1
    ##   Cod       1      1
    ##   Saithe    1      1

For the North Sea this is not the case and so the model would be
improved by also including an interaction matrix which describes the
spatial overlap between species.

An example interaction matrix for the North Sea has been included in
`mizer` as a .csv file. The location of the file can be found by
running:

``` r

inter_location <- system.file("extdata", "NS_interaction.csv",
                              package = "mizer")
```

Take a look at it in a spreadsheet if you want. As mentioned above, to
read this file into R we can make use of the
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) function.
However, this time we want the first column of the .csv file to be the
row names. We therefore use an additional argument to the
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) function:
`row.names`.

``` r

inter <- read.csv(inter_location, row.names = 1)
inter
```

    ##              Sprat    Sandeel     N.pout    Herring        Dab    Whiting
    ## Sprat   0.72912919 0.03408440 0.06354517 0.27416982 0.36241552 0.26525924
    ## Sandeel 0.03408440 0.68119882 0.04892432 0.05888214 0.09736663 0.07510011
    ## N.pout  0.06354517 0.04892432 0.79660429 0.29755069 0.09088798 0.29989886
    ## Herring 0.27416982 0.05888214 0.29755069 0.65890104 0.28963957 0.37373927
    ## Dab     0.36241552 0.09736663 0.09088798 0.28963957 0.80817768 0.33389727
    ## Whiting 0.26525924 0.07510011 0.29989886 0.37373927 0.33389727 0.70928230
    ## Sole    0.29795558 0.06020860 0.01679020 0.20014139 0.38047464 0.19227455
    ## Gurnard 0.17515576 0.05992649 0.30624141 0.27510627 0.22041200 0.37109904
    ## Plaice  0.37065975 0.07801855 0.07855818 0.27791867 0.56492206 0.29503807
    ## Haddock 0.08135547 0.09395730 0.54917554 0.34835469 0.13168065 0.39164787
    ## Cod     0.33757969 0.09943453 0.32502256 0.40477930 0.41647801 0.44060879
    ## Saithe  0.01681321 0.01609022 0.29498937 0.12620591 0.03138197 0.10228168
    ##               Sole    Gurnard     Plaice    Haddock        Cod     Saithe
    ## Sprat   0.29795558 0.17515576 0.37065975 0.08135547 0.33757969 0.01681321
    ## Sandeel 0.06020860 0.05992649 0.07801855 0.09395730 0.09943453 0.01609022
    ## N.pout  0.01679020 0.30624141 0.07855818 0.54917554 0.32502256 0.29498937
    ## Herring 0.20014139 0.27510627 0.27791867 0.34835469 0.40477930 0.12620591
    ## Dab     0.38047464 0.22041200 0.56492206 0.13168065 0.41647801 0.03138197
    ## Whiting 0.19227455 0.37109904 0.29503807 0.39164787 0.44060879 0.10228168
    ## Sole    0.71558049 0.10677895 0.39137317 0.03447799 0.25761229 0.01242055
    ## Gurnard 0.10677895 0.88010500 0.16492120 0.35735444 0.35183282 0.12351994
    ## Plaice  0.39137317 0.16492120 0.71922391 0.11248513 0.35043671 0.03294939
    ## Haddock 0.03447799 0.35735444 0.11248513 0.85830725 0.39577341 0.26167470
    ## Cod     0.25761229 0.35183282 0.35043671 0.39577341 0.78654705 0.20894496
    ## Saithe  0.01242055 0.12351994 0.03294939 0.26167470 0.20894496 0.66383553

We can set the interaction matrix in our existing MizerParams object
`params` with the
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md)
function:

``` r

params <- setInteraction(params, interaction = inter)
```

Alternatively, instead of changing the interaction matrix in the
existing MizerParams object, we could have created a new object from
scratch with our interaction matrix by passing it to
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md):

``` r

params_new <- newMultispeciesParams(species_params, interaction = inter)
```

    ## Because you have n != p, the default value for `h` is not very good.
    ## Because the age at maturity is not known, I need to fall back to using
    ## von Bertalanffy parameters, where available, and this is not reliable.
    ## No ks column so calculating from critical feeding level.
    ## Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.
    ## Using f0, h, lambda, kappa and the predation kernel to calculate gamma.

Note that the first argument must be the species parameters data.frame.
The remaining arguments can be in any order but should be named. We are
using the default values for all other parameters.

We now have all we need to start running projections. Before we get to
that though, we’ll take a quick look at how different fishing gears can
be set up.

### Setting different gears

In the above example, each species is caught by the same gear (named
“knife_edge_gear”). This is the default when no gear information is
provided.

``` r

gear_params(params)
```

    ##                          species            gear   sel_func knife_edge_size
    ## Sprat, knife_edge_gear     Sprat knife_edge_gear knife_edge              13
    ## Sandeel, knife_edge_gear Sandeel knife_edge_gear knife_edge               4
    ## N.pout, knife_edge_gear   N.pout knife_edge_gear knife_edge              23
    ## Herring, knife_edge_gear Herring knife_edge_gear knife_edge              99
    ## Dab, knife_edge_gear         Dab knife_edge_gear knife_edge              21
    ## Whiting, knife_edge_gear Whiting knife_edge_gear knife_edge              75
    ## Sole, knife_edge_gear       Sole knife_edge_gear knife_edge              78
    ## Gurnard, knife_edge_gear Gurnard knife_edge_gear knife_edge              39
    ## Plaice, knife_edge_gear   Plaice knife_edge_gear knife_edge             105
    ## Haddock, knife_edge_gear Haddock knife_edge_gear knife_edge             165
    ## Cod, knife_edge_gear         Cod knife_edge_gear knife_edge            1606
    ## Saithe, knife_edge_gear   Saithe knife_edge_gear knife_edge            1076
    ##                          catchability
    ## Sprat, knife_edge_gear              1
    ## Sandeel, knife_edge_gear            1
    ## N.pout, knife_edge_gear             1
    ## Herring, knife_edge_gear            1
    ## Dab, knife_edge_gear                1
    ## Whiting, knife_edge_gear            1
    ## Sole, knife_edge_gear               1
    ## Gurnard, knife_edge_gear            1
    ## Plaice, knife_edge_gear             1
    ## Haddock, knife_edge_gear            1
    ## Cod, knife_edge_gear                1
    ## Saithe, knife_edge_gear             1

Here, we look at an example where we set up four different gears:
Industrial, Pelagic, Beam and Otter trawl, that catch different
combinations of species. We can achieve that by only changing the `gear`
column in the `gear_params` data frame.

``` r

gear_params(params)$gear <- c("Industrial", "Industrial", "Industrial",
                              "Pelagic", "Beam", "Otter",
                              "Beam", "Otter", "Beam",
                              "Otter", "Otter", "Otter")
```

You can see the result by calling
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md) on
the `params` object.

``` r

summary(params)
```

    ## An object of class "MizerParams" 
    ## mizer version: 3.1.0
    ## Created: 2026-06-25 21:02:52
    ## Modified: 2026-06-25 21:02:52
    ## Consumer size spectrum:
    ##  minimum size:   0.001
    ##  maximum size:   60066
    ##  no. size bins:  100
    ## Resource size spectrum:
    ##  minimum size:   9.20914e-13
    ##  maximum size:   8.48399
    ##  no. size bins:  166 (215 size bins in total)
    ## Species details:
    ##    species w_inf w_mat w_min   beta sigma
    ## 1    Sprat    33    13 0.001  51076   0.8
    ## 2  Sandeel    36     4 0.001 398849   1.9
    ## 3   N.pout   100    23 0.001     22   1.5
    ## 4  Herring   334    99 0.001 280540   3.2
    ## 5      Dab   324    21 0.001    191   1.9
    ## 6  Whiting  1192    75 0.001     22   1.5
    ## 7     Sole   866    78 0.001    381   1.9
    ## 8  Gurnard   668    39 0.001    283   1.8
    ## 9   Plaice  2976   105 0.001    113   1.6
    ## 10 Haddock  3485   165 0.001    558   2.1
    ## 11     Cod 40044  1606 0.001     66   1.3
    ## 12  Saithe 16856  1076 0.001     40   1.1
    ## 
    ## Fishing gear details:
    ## Gear          Effort  Target species 
    ##  ----------------------------------
    ## Industrial     0.00   Sprat, Sandeel, N.pout 
    ## Pelagic        0.00   Herring 
    ## Beam           0.00   Dab, Sole, Plaice 
    ## Otter          0.00   Whiting, Gurnard, Haddock, Cod, Saithe

In this example the same gear now catches multiple stocks. For example,
the *Industrial* gear catches Sprat, Sandeel and Norway Pout. Why would
we want to set up the gears like this? In the next section on [running a
multispecies
model](https://sizespectrum.org/mizer/articles/running_a_simulation.md)
we will see that to project the model through time you can specify the
fishing effort for each gear through time. By setting the gears up in
this way you can run different management scenarios of changing the
efforts of the fishing gears rather than on individual species. It also
means that after a simulation has been run you can examine the catches
by gear.

### Setting to steady state

Once the `MizerParams` object has been properly set up, it may be the
case that one wishes put the system in steady state. Sometimes this can
be done simply by running the model using
[`project()`](https://sizespectrum.org/mizer/reference/project.md) until
it reaches steady state. However, this method is not guaranteed to work,
and there is a function called
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) that is
more reliable. The function
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) must be
supplied with a MizerParams object. It takes that MizerParams object,
looks at the initial system state, computes the levels of reproduction
of the different species, hold them fixed, and evolves the system until
a steady state is reached (or more precisely, until the amount that the
population abundances change during a time-step is below some small
tolerance level). After this, the reproductive efficiency of each
species is altered so that when the reproduction dynamics are turned
back on (i.e., when we stop holding recruitment levels fixed), the
values of the reproduction levels which we held the system fixed at will
be realized. The steady function is not sure to converge, and the way it
re-tunes the reproductive efficiency values may not be realistic, but
the idea is to alter the other parameters in the system until
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) does
arrive at a steady state with sensible reproductive efficiency values.

Now that we know how to create a multispecies model we shall discuss how
to [run a multispecies
model.](https://sizespectrum.org/mizer/articles/running_a_simulation.md)
