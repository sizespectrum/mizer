# Set fishing parameters

Set fishing parameters

## Usage

``` r
setFishing(
  params,
  selectivity = NULL,
  catchability = NULL,
  reset = FALSE,
  initial_effort = NULL,
  ...
)

getCatchability(params)

catchability(params)

catchability(params) <- value

getSelectivity(params)

selectivity(params)

selectivity(params) <- value

getInitialEffort(params)
```

## Arguments

- params:

  A MizerParams object

- selectivity:

  Optional. An array (gear x species x size) that holds the selectivity
  of each gear for species and size, \\S\_{g,i,w}\\.

- catchability:

  Optional. An array (gear x species) that holds the catchability of
  each species by each gear, \\Q\_{g,i}\\.

- reset:

  If set to TRUE, then both `catchability` and `selectivity` will be
  reset to the values calculated from the gear parameters, even if it
  was previously overwritten with a custom value. If set to FALSE
  (default) then a recalculation from the gear parameters will take
  place only if no custom value has been set.

- initial_effort:

  Optional. A number or a named numeric vector specifying the fishing
  effort. If a number, the same effort is used for all gears. If a
  vector, must be named by gear.

- ...:

  Unused

- value:

  The array to assign

## Value

`setFishing()`: A MizerParams object with updated fishing parameters.

`getCatchability()` or equivalently `catchability()`: An array (gear x
species) that holds the catchability of each species by each gear,
\\Q\_{g,i}\\. The names of the dimensions are "gear, "sp".

`getSelectivity()` or equivalently `selectivity()`: An array (gear x
species x size) that holds the selectivity of each gear for species and
size, \\S\_{g,i,w}\\. The names of the dimensions are "gear, "sp", "w".

`getInitialEffort()` or equivalently
[`initial_effort()`](https://sizespectrum.org/mizer/reference/initial_effort.md):
A named vector with the initial fishing effort for each gear.

## Setting fishing

**Gears**

In `mizer`, fishing mortality is imposed on species by fishing gears.
The total per-capita fishing mortality (1/year) is obtained by summing
over the mortality from all gears, \$\$\mu\_{f.i}(w) = \sum_g
F\_{g,i}(w),\$\$ where the fishing mortality \\F\_{g,i}(w)\\ imposed by
gear \\g\\ on species \\i\\ at size \\w\\ is calculated as:
\$\$F\_{g,i}(w) = S\_{g,i}(w) Q\_{g,i} E\_{g},\$\$ where \\S\\ is the
selectivity by species, gear and size, \\Q\\ is the catchability by
species and gear and \\E\\ is the fishing effort by gear.

**Selectivity**

The selectivity at size of each gear for each species is saved as a
three dimensional array (gear x species x size). Each entry has a range
between 0 (that gear is not selecting that species at that size) to 1
(that gear is selecting all individuals of that species of that size).
This three dimensional array can be specified explicitly via the
`selectivity` argument, but usually mizer calculates it from the
`gear_params` slot of the MizerParams object.

To allow the calculation of the `selectivity` array, the `gear_params`
slot must be a data frame with one row for each gear-species
combination. So if for example a gear can select three species, then
that gear contributes three rows to the `gear_params` data frame, one
for each species it can select. The data frame must have columns `gear`,
holding the name of the gear, `species`, holding the name of the
species, and `sel_func`, holding the name of the function that
calculates the selectivity curve. Some selectivity functions are
included in the package:
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md),
[`sigmoid_length()`](https://sizespectrum.org/mizer/reference/sigmoid_length.md),
[`double_sigmoid_length()`](https://sizespectrum.org/mizer/reference/double_sigmoid_length.md),
and
[`sigmoid_weight()`](https://sizespectrum.org/mizer/reference/sigmoid_weight.md).
Users are able to write their own size-based selectivity function. The
first argument to the function must be `w` and the function must return
a vector of the selectivity (between 0 and 1) at size.

Each selectivity function may have parameters. Values for these
parameters must be included as columns in the gear parameters
data.frame. The names of the columns must exactly match the names of the
corresponding arguments of the selectivity function. For example, the
default selectivity function is
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md)
that a has sudden change of selectivity from 0 to 1 at a certain size.
In its help page you can see that the
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md)
function has arguments `w` and `knife_edge_size`. The first argument,
`w`, is size (the function calculates selectivity at size). All
selectivity functions must have `w` as the first argument. The values
for the other arguments must be found in the gear parameters data.frame.
So for the
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md)
function there should be a `knife_edge_size` column. Because
[`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md)
is the default selectivity function, the `knife_edge_size` argument has
a default value = `w_mat`.

The most commonly-used selectivity function is
[`sigmoid_length()`](https://sizespectrum.org/mizer/reference/sigmoid_length.md).
It has a smooth transition from 0 to 1 at a certain size. The
[`sigmoid_length()`](https://sizespectrum.org/mizer/reference/sigmoid_length.md)
function has the two parameters `l50` and `l25` that are the lengths in
cm at which 50% or 25% of the fish are selected by the gear. If you
choose this selectivity function then the `l50` and `l25` columns must
be included in the gear parameters data.frame.

In case each species is only selected by one gear, the columns of the
`gear_params` data frame can alternatively be provided as columns of the
`species_params` data frame, if this is more convenient for the user to
set up. Mizer will then copy these columns over to create the
`gear_params` data frame when it creates the MizerParams object. However
changing these columns in the species parameter data frame later will
not update the `gear_params` data frame.

**Catchability**

Catchability is used as an additional factor to make the link between
gear selectivity, fishing effort and fishing mortality. For example, it
can be set so that an effort of 1 gives a desired fishing mortality. In
this way effort can then be specified relative to a 'base effort', e.g.
the effort in a particular year.

Catchability is stored as a two dimensional array (gear x species). This
can either be provided explicitly via the `catchability` argument, or
the information can be provided via a `catchability` column in the
`gear_params` data frame.

In the case where each species is selected by only a single gear, the
`catchability` column can also be provided in the `species_params` data
frame. Mizer will then copy this over to the `gear_params` data frame
when the MizerParams object is created.

**Effort**

The initial fishing effort is stored in the `MizerParams` object. If it
is not supplied, it is set to zero. The initial effort can be overruled
when the simulation is run with
[`project()`](https://sizespectrum.org/mizer/reference/project.md),
where it is also possible to specify an effort that varies through time.

## See also

[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)

Other functions for setting parameters:
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
[`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
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
# Halve the initial fishing effort for all gears
params <- setFishing(NS_params, initial_effort = 0.5)
getInitialEffort(params)
#> Industrial    Pelagic       Beam      Otter 
#>        0.5        0.5        0.5        0.5 
str(getCatchability(NS_params))
#>  num [1:4, 1:12] 1 0 0 0 1 0 0 0 1 0 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ gear: chr [1:4] "Industrial" "Pelagic" "Beam" "Otter"
#>   ..$ sp  : chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
str(getSelectivity(NS_params))
#>  num [1:4, 1:12, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ gear: chr [1:4] "Industrial" "Pelagic" "Beam" "Otter"
#>   ..$ sp  : chr [1:12] "Sprat" "Sandeel" "N.pout" "Herring" ...
#>   ..$ w   : chr [1:100] "0.001" "0.00119" "0.00142" "0.0017" ...
str(getInitialEffort(NS_params))
#>  Named num [1:4] 0 1 0.5 0.5
#>  - attr(*, "names")= chr [1:4] "Industrial" "Pelagic" "Beam" "Otter"
```
