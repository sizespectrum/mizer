# Set up parameters for a general multispecies model

Sets up a multi-species size spectrum model by filling all slots in the
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
object based on user-provided or default parameters. There is a long
list of arguments, but almost all of them have sensible default values.
The only required argument is the `species_params` data frame. All
arguments are described in more details in the sections below the list.

## Usage

``` r
newMultispeciesParams(
  species_params,
  interaction = NULL,
  no_w = 100,
  min_w = 0.001,
  max_w = NA,
  min_w_pp = NA,
  pred_kernel = NULL,
  search_vol = NULL,
  intake_max = NULL,
  metab = NULL,
  p = 0.7,
  ext_mort = NULL,
  z0pre = 0.6,
  z0exp = n - 1,
  ext_encounter = NULL,
  maturity = NULL,
  repro_prop = NULL,
  RDD = "BevertonHoltRDD",
  kappa = 1e+11,
  n = 2/3,
  resource_rate = 10,
  resource_capacity = kappa,
  lambda = 2.05,
  w_pp_cutoff = 10,
  resource_dynamics = "resource_semichemostat",
  gear_params = NULL,
  selectivity = NULL,
  catchability = NULL,
  initial_effort = NULL,
  second_order_w = FALSE,
  info_level = 3,
  z0 = deprecated(),
  r_pp = deprecated()
)
```

## Arguments

- species_params:

  A data frame of species-specific parameter values.

- interaction:

  Optional interaction matrix of the species (predator species x prey
  species). By default all entries are 1. See "Setting interaction
  matrix" section below.

- no_w:

  The number of size bins in the consumer spectrum.

- min_w:

  Sets the size of the eggs of all species for which this is not given
  in the `w_min` column of the `species_params` dataframe.

- max_w:

  The largest size of the consumer spectrum. By default this is set to
  the largest `w_max` specified in the `species_params` data frame.

- min_w_pp:

  The smallest size of the resource spectrum. By default this is set to
  the smallest value at which any of the consumers can feed.

- pred_kernel:

  Optional. An array (species x predator size x prey size) that holds
  the predation coefficient of each predator at size on each prey size.
  If not supplied, a default is set as described in section "Setting
  predation kernel".

- search_vol:

  Optional. An array (species x size) holding the search volume for each
  species at size. If not supplied, a default is set as described in the
  section "Setting search volume".

- intake_max:

  Optional. An array (species x size) holding the maximum intake rate
  for each species at size. If not supplied, a default is set as
  described in the section "Setting maximum intake rate".

- metab:

  Optional. An array (species x size) holding the metabolic rate for
  each species at size. If not supplied, a default is set as described
  in the section "Setting metabolic rate".

- p:

  The allometric metabolic exponent. This is only used if `metab` is not
  given explicitly and if the exponent is not specified in a `p` column
  in the `species_params`.

- ext_mort:

  Optional. An array (species x size) holding the external mortality
  rate. If not supplied, a default is set as described in the section
  "Setting external mortality rate".

- z0pre:

  If `z0`, the mortality from other sources, is not a column in the
  species data frame, it is calculated as z0pre \* w_inf ^ z0exp.
  Default value is 0.6.

- z0exp:

  If `z0`, the mortality from other sources, is not a column in the
  species data frame, it is calculated as `z0pre * w_inf ^ z0exp`.
  Default value is `n-1`.

- ext_encounter:

  Optional. An array (species x size) holding the external encounter
  rate. If not supplied, a default is calculated from the `E_ext` and
  `n` species parameters as described in the section "Setting external
  encounter rate".

- maturity:

  Optional. An array (species x size) that holds the proportion of
  individuals of each species at size that are mature. If not supplied,
  a default is set as described in the section "Setting reproduction".

- repro_prop:

  Optional. An array (species x size) that holds the proportion of the
  energy available for growth and reproduction that a mature individual
  allocates to reproduction for each species at size. If not supplied, a
  default is set as described in the section "Setting reproduction".

- RDD:

  The name of the function calculating the density-dependent
  reproduction rate from the density-independent rate. Defaults to
  "[`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md)".

- kappa:

  The coefficient of the initial resource abundance power-law.

- n:

  The allometric growth exponent. This can be overruled for individual
  species by including a `n` column in the `species_params`.

- resource_rate:

  Optional. A vector of per-capita resource birth rate for each size
  class or a single number giving the coefficient in the power-law for
  this rate, see "Setting resource dynamics" below. Must be strictly
  positive.

- resource_capacity:

  Optional. Vector of resource intrinsic carrying capacities or
  coefficient in the power-law for the capacity, see "Setting resource
  dynamics" below. The resource capacity must not be smaller than the
  resource abundance.

- lambda:

  Used to set power-law exponent for resource capacity if the
  `resource_capacity` argument is given as a single number.

- w_pp_cutoff:

  The upper cut off size of the resource spectrum power law used when
  `resource_capacity` is given as a single number. When changing
  `w_pp_cutoff` without providing `resource_capacity`, the cutoff can
  only be decreased. In that case, both the carrying capacity and the
  initial resource abundance will be cut off at the new value. To
  increase the cutoff, you must also provide the `resource_capacity` for
  the extended range.

- resource_dynamics:

  Optional. Name of the function that determines the resource dynamics
  by calculating the resource spectrum at the next time step from the
  current state.

- gear_params:

  A data frame with gear-specific parameter values.

- selectivity:

  Optional. An array (gear x species x size) that holds the selectivity
  of each gear for species and size, \\S\_{g,i,w}\\.

- catchability:

  Optional. An array (gear x species) that holds the catchability of
  each species by each gear, \\Q\_{g,i}\\.

- initial_effort:

  Optional. A number or a named numeric vector specifying the fishing
  effort. If a number, the same effort is used for all gears. If a
  vector, must be named by gear.

- second_order_w:

  **\[experimental\]** Selects the second-order numerical scheme for the
  new model. Accepts the same values as the
  [`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
  setter: a single logical (`TRUE` switches on both second-order flux
  and bin-averaging), a single flux scheme name (`"upwind"`,
  `"van_leer"` or `"centred"`), or a named vector with entries `flux`
  and/or `bin_average`. The `bin_average` choice is applied *before* the
  resource and abundance power laws are constructed, so they are built
  bin-averaged from the start (unlike setting
  [`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
  on an existing object). The `flux` scheme governs time projection
  only, so the robust first-order upwind scheme is used for the
  construction-time steady-state solve and the chosen scheme is then
  activated for the returned model. Defaults to `FALSE` (the first-order
  behaviour of previous mizer).

- info_level:

  Controls the amount of information messages that are shown when the
  function sets default values for parameters. Higher levels lead to
  more messages.

- z0:

  **\[deprecated\]** Use `ext_mort` instead. Not to be confused with the
  species_parameter `z0`.

- r_pp:

  **\[deprecated\]**. Use `resource_rate` argument instead.

## Value

An object of type
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)

## Species parameters

The only essential argument is a data frame that contains the species
parameters. The data frame is arranged species by parameter, so each
column of the parameter data frame is a parameter and each row has the
values of the parameters for one of the species in the model.

There are two essential columns that must be included in the species
parameter data.frame and that do not have default values: the `species`
column that should hold strings with the names of the species and the
`w_inf` column with the von Bertalanffy asymptotic sizes of the species
in grams. (You could alternatively specify the corresponding length in
cm in an `l_inf` column.) The computational upper size boundary `w_max`
is not essential; if it is missing it defaults to `1.5 * w_inf`. For
backwards compatibility, if `w_inf` is missing it is taken from the
`w_repro_max` or `w_max` column instead.

The `species_params dataframe` also needs to contain the parameters
needed by any predation kernel function (size selectivity function).
This will be mentioned in the appropriate sections below.

For all other species parameters, mizer will calculate default values if
they are not included in the species parameter data frame. They will be
automatically added when the `MizerParams` object is created. For these
parameters you can also specify values for only some species and leave
the other entries as NA and the missing values will be set to the
defaults. So the `species_params` data frame saved in the returned
MizerParams object will differ from the one you supply because it will
have the missing species parameters filled in with default values.

If you are not happy with any of the species parameter values used you
can always change them later with
[`species_params<-()`](https://sizespectrum.org/mizer/reference/species_params.md).

All the parameters will be mentioned in the following sections.

## Setting initial values

The initial values for the species number densities are set using the
function
[`get_initial_n()`](https://sizespectrum.org/mizer/reference/get_initial_n.md).
These are quite arbitrary and not very close to the steady state
abundances. We intend to improve this in the future.

The initial resource number density \\N_R(w)\\ is set to a power law
with coefficient `kappa` (\\\kappa\\) and exponent `-lambda`
(\\-\lambda\\): \$\$N_R(w) = \kappa\\ w^{-\lambda}\$\$ for all \\w\\
less than `w_pp_cutoff` and zero for sizes at or above `w_pp_cutoff`.

## Size grid

A size grid is created so that the log-sizes are equally spaced. The
spacing is chosen so that there will be `no_w` fish size bins, with the
smallest starting at `min_w` and the largest starting at `max_w`. For
the resource spectrum there is a larger set of bins containing
additional bins below `min_w`, with the same log size. The number of
extra bins is such that `min_w_pp` comes to lie within the smallest bin.

## Units in mizer

Mizer uses grams to measure weight, centimetres to measure lengths, and
years to measure time.

Mizer is agnostic about whether abundances are given as

1.  numbers per area,

2.  numbers per volume or

3.  total numbers for the entire study area.

You should make the choice most convenient for your application and then
stick with it. If you make choice 1 or 2 you will also have to choose a
unit for area or volume. Your choice will then determine the units for
some of the parameters. This will be mentioned when the parameters are
discussed in the sections below.

Your choice will also affect the units of the quantities you may want to
calculate with the model. For example, the yield will be in
grams/year/m^2 in case 1 if you choose m^2 as your measure of area, in
grams/year/m^3 in case 2 if you choose m^3 as your unit of volume, or
simply grams/year in case 3. The same comment applies for other
measures, like total biomass, which will be grams/area in case 1,
grams/volume in case 2 or simply grams in case 3. When mizer puts units
on axes in plots, it will choose the units appropriate for case 3. So
for example in
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md)
it gives the unit as grams.

You can convert between these choices. For example, if you use case 1,
you need to multiply with the area of the ecosystem to get the total
quantity. If you work with case 2, you need to multiply by both area and
the thickness of the productive layer. In that respect, case 2 is a bit
cumbersome. The function
[`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.md)
is useful to change the units you are using.

## Setting interaction matrix

You do not need to specify an interaction matrix. If you do not, then
the predator-prey interactions are purely determined by the size of
predator and prey and totally independent of the species of predator and
prey.

The interaction matrix \\\theta\_{ij}\\ modifies the interaction of each
pair of species in the model. This can be used for example to allow for
different spatial overlap among the species. The values in the
interaction matrix are used to scale the encountered food and predation
mortality (see on the website [the section on predator-prey encounter
rate](https://sizespectrum.org/mizer/articles/model_description.html#sec:pref)
and on [predation
mortality](https://sizespectrum.org/mizer/articles/model_description.html#mortality)).
The first index refers to the predator species and the second to the
prey species.

The interaction matrix is used when calculating the food encounter rate
in
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
and the predation mortality rate in
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md).
Its entries are dimensionless numbers. If all the values in the
interaction matrix are equal then predator-prey interactions are
determined entirely by size-preference.

This function checks that the supplied interaction matrix is valid and
then stores it in the `interaction` slot of the `params` object.

The order of the columns and rows of the `interaction` argument should
be the same as the order in the species params data frame in the
`params` object. If you supply a named array then the function will
check the order and message if it is different before ignoring the
supplied dimnames. If you supply only column names then these are also
used as the row names. One way of creating your own interaction matrix
is to enter the data using a spreadsheet program and saving it as a .csv
file. The data can then be read into R using the command
[`read.csv()`](https://rdrr.io/r/utils/read.table.html).

The interaction of the species with the resource are set via a column
`interaction_resource` in the `species_params` data frame. By default
this column is set to all 1s.

## Setting predation kernel

**Kernel dependent on predator to prey size ratio**

If the `pred_kernel` argument is not supplied, then this function sets a
predation kernel that depends only on the ratio of predator mass to prey
mass, not on the two masses independently. The shape of that kernel is
then determined by the `pred_kernel_type` column in species_params.

The default for `pred_kernel_type` is "lognormal". This will call the
function
[`lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/lognormal_pred_kernel.md)
to calculate the predation kernel. An alternative pred_kernel type is
"box", implemented by the function
[`box_pred_kernel()`](https://sizespectrum.org/mizer/reference/box_pred_kernel.md),
and "power_law", implemented by the function
[`power_law_pred_kernel()`](https://sizespectrum.org/mizer/reference/power_law_pred_kernel.md).
These functions require certain species parameters in the species_params
data frame. For the lognormal kernel these are `beta` and `sigma`, for
the box kernel they are `ppmr_min` and `ppmr_max`. They are explained in
the help pages for the kernel functions. Except for `beta` and `sigma`,
no defaults are set for these parameters. If they are missing from the
species_params data frame then mizer will issue an error message.

You can use any other string for `pred_kernel_type`. If for example you
choose "my" then you need to define a function `my_pred_kernel` that you
can model on the existing functions like
[`lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/lognormal_pred_kernel.md).

When using a kernel that depends on the predator/prey size ratio only,
mizer does not need to store the entire three dimensional array in the
MizerParams object. Such an array can be very big when there is a large
number of size bins. Instead, mizer only needs to store two
two-dimensional arrays that hold Fourier transforms of the feeding
kernel function that allow the encounter rate and the predation rate to
be calculated very efficiently. However, if you need the full
three-dimensional array you can calculate it with the
[`getPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
function.

**Kernel dependent on both predator and prey size**

If you want to work with a feeding kernel that depends on predator mass
and prey mass independently, you can specify the full feeding kernel as
a three-dimensional array (predator species x predator size x prey
size).

You should use this option only if a kernel dependent only on the
predator/prey mass ratio is not appropriate. Using a kernel dependent on
predator/prey mass ratio only allows mizer to use fast Fourier transform
methods to significantly reduce the running time of simulations.

The order of the predator species in `pred_kernel` should be the same as
the order in the species params dataframe in the `params` object. If you
supply a named array then the function will check the order and warn if
it is different.

## Setting search volume

The search volume \\\gamma_i(w)\\ of an individual of species \\i\\ and
weight \\w\\ multiplies the predation kernel when calculating the
encounter rate in
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
and the predation rate in
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md).

The name "search volume" is a bit misleading, because \\\gamma_i(w)\\
does not have units of volume. It is simply a parameter that determines
the rate of predation. Its units depend on your choice, see section
"Units in mizer". If you have chosen to work with total abundances, then
it is a rate with units 1/year. If you have chosen to work with
abundances per m^2 then it has units of m^2/year. If you have chosen to
work with abundances per m^3 then it has units of m^3/year.

If the `search_vol` argument is not supplied, then the search volume is
set to \$\$\gamma_i(w) = \gamma_i w^q_i.\$\$ The values of \\\gamma_i\\
(the search volume at 1g) and \\q_i\\ (the allometric exponent of the
search volume) are taken from the `gamma` and `q` columns in the species
parameter dataframe. If the `gamma` column is not supplied in the
species parameter dataframe, a default is calculated by the
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md)
function. If the `q` column is not supplied, a default of
`lambda - 2 + n` is used. Note that only for predators of size \\w = 1\\
gram is the value of the species parameter \\\gamma_i\\ the same as the
value of the search volume \\\gamma_i(w)\\.

If the `search_vol` slot has a comment and `reset = FALSE`, then a
recalculation from the species parameters is suppressed and a message is
issued if the recalculated values would differ from the stored ones.

## Setting maximum intake rate

The maximum intake rate \\h_i(w)\\ of an individual of species \\i\\ and
weight \\w\\ determines the feeding level, calculated with
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md).
It is measured in grams/year.

If the `intake_max` argument is not supplied, then the maximum intake
rate is set to \$\$h_i(w) = h_i w^{n_i}.\$\$ The values of \\h_i\\ (the
maximum intake rate of an individual of size 1 gram) and \\n_i\\ (the
allometric exponent for the intake rate) are taken from the `h` and `n`
columns in the species parameter dataframe. If the `h` column is not
supplied in the species parameter dataframe, it is calculated by the
[`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md)
function. If the `n` column is not supplied, a default of \\n_i = 3/4\\
is used.

If \\h_i\\ is set to `Inf`, fish of species i will consume all
encountered food.

If the `intake_max` slot has a comment and `reset = FALSE`, then a
recalculation from the species parameters is suppressed and a message is
issued if the recalculated values would differ from the stored ones.

## Setting metabolic rate

The metabolic rate is subtracted from the energy income rate to
calculate the rate at which energy is available for growth and
reproduction, see
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md).
It is measured in grams/year.

If the `metab` argument is not supplied, then for each species the
metabolic rate \\k(w)\\ for an individual of size \\w\\ is set to
\$\$k(w) = k_s w^p + k w,\$\$ where \\k_s w^p\\ represents the rate of
standard metabolism and \\k w\\ is the rate at which energy is expended
on activity and movement. The values of \\k_s\\, \\p\\ and \\k\\ are
taken from the `ks`, `p` and `k` columns in the species parameter
dataframe. If any of these parameters are not supplied, the defaults are
\\k = 0\\, \\p = 3/4\\ and \$\$k_s = f_c h \alpha w\_{mat}^{n-p},\$\$
where \\f_c\\ is the critical feeding level taken from the `fc` column
in the species parameter data frame. If the critical feeding level is
not specified, a default of \\f_c = 0.2\\ is used.

If the `metab` slot has a comment and `reset = FALSE`, then a
recalculation from the species parameters is suppressed and a message is
issued if the recalculated values would differ from the stored ones.

## Setting external mortality rate

The external mortality is all the mortality that is not due to fishing
or predation by predators included in the model. The external mortality
could be due to predation by predators that are not explicitly included
in the model (e.g. mammals or seabirds) or due to other causes like
illness. It is a rate with units 1/year.

The `ext_mort` argument allows you to specify an external mortality rate
that depends on species and body size. You can see an example of this in
the Examples section of the help page for
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md).

If the `ext_mort` argument is not supplied, then the external mortality
is taken from the species parameters as \$\$\mu\_{ext.i}(w) = z\_{0.i} +
z\_{ext.i} w^{d_i}.\$\$ The value of the constant \\z_0\\ for each
species is taken from the `z0` column of the species parameter data
frame, if that column exists. Otherwise it is calculated as \$\$z\_{0.i}
= {\tt z0pre}\_i\\ w\_{inf}^{\tt z0exp}.\$\$ Missing values of `z_ext`
are set to 0 and missing values of `d` are set to `n - 1`.

By default the power law is evaluated at the left bin edges \\w_j\\
(point sampling). If the `bin_average` entry of the `second_order_w`
slot is `TRUE` (see
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)),
then the \\z\_{ext} w^d\\ term is instead replaced by its exact average
over each bin \\\[w_j, w\_{j+1}\]\\, \$\$\frac{z\_{ext}}{\Delta
w_j}\int\_{w_j}^{w\_{j+1}} w^d\\ dw = z\_{ext}\\\frac{w\_{j+1}^{d+1} -
w_j^{d+1}}{(d+1)\\\Delta w_j},\$\$ (with the limiting form
\\z\_{ext}\ln(w\_{j+1}/w_j)/\Delta w_j\\ when \\d = -1\\). This is the
consistent choice in the finite-volume scheme, where the external
mortality multiplies the bin-averaged abundance. The bin-averaging is
applied only to the auto-calculated power-law default; a user-supplied
`ext_mort` array is left untouched.

## Setting external encounter rate

The external encounter rate is the rate at which a predator encounters
food that is not explicitly modelled. It is a rate with units mass/year.

The `ext_encounter` argument allows you to specify an external encounter
rate that depends on species and body size. You can see an example of
this in the Examples section of the help page for
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md).

If the `ext_encounter` argument is not supplied, then the external
encounter rate is calculated as a power law: \$\$E\_{ext.i}(w) =
E\_{ext.i}\\ w^{n_i}.\$\$ The coefficient \\E\_{ext.i}\\ is taken from
the `E_ext` column of the species parameter data frame, which defaults
to 0. The exponent \\n_i\\ is taken from the `n` column of the species
parameter data frame.

If the `ext_encounter` slot has a comment and `reset = FALSE`, then a
recalculation from the species parameters is suppressed and a message is
issued if the recalculated values would differ from the stored ones.

## Setting external diffusion rate

The external diffusion rate allows you to impose additional diffusion
beyond the predation-driven diffusion that can be internally modelled by
mizer.

The `ext_diffusion` argument allows you to specify a diffusion rate that
depends on species and body size.

If the `ext_diffusion` argument is not supplied, then the external
diffusion rate is calculated as a power law: \$\$D\_{ext.i}(w) =
D\_{ext.i}\\ w^{n_i+1}.\$\$ The coefficient \\D\_{ext.i}\\ is taken from
the `D_ext` column of the species parameter data frame, which defaults
to 0. The exponent \\n_i + 1\\ uses the `n` column of the species
parameter data frame.

If the `ext_diffusion` slot has a comment and `reset = FALSE`, then a
recalculation from the species parameters is suppressed and a message is
issued if the recalculated values would differ from the stored ones.

## Setting reproduction

For each species and at each size, the proportion \\\psi\\ of the
available energy that is invested into reproduction is the product of
two factors: the proportion `maturity` of individuals that are mature
and the proportion `repro_prop` of the energy available to a mature
individual that is invested into reproduction. There is a size
`w_repro_max` at which a typical mature individual invests all of its
available energy into reproduction. This is **not** a hard ceiling on
size: not all individuals are mature at `w_repro_max`, and diffusion in
the growth process allows some individuals to grow beyond it, so fish
larger than `w_repro_max` can exist. If you have not specified the
`w_repro_max` column in the species parameter data frame, then the von
Bertalanffy asymptotic size `w_inf` is used instead.

### Maturity ogive

If the the proportion of individuals that are mature is not supplied via
the `maturity` argument, then it is set to a sigmoidal maturity ogive
that changes from 0 to 1 at around the maturity size: \$\${\tt
maturity}(w) =
\left\[1+\left(\frac{w}{w\_{mat}}\right)^{-U}\right\]^{-1}.\$\$ (To
avoid clutter, we are not showing the species index in the equations,
although each species has its own maturity ogive.) The maturity weights
are taken from the `w_mat` column of the species_params data frame. Any
missing maturity weights are set to 1/4 of the asymptotic size in the
`w_inf` column.

The exponent \\U\\ determines the steepness of the maturity ogive. By
default it is chosen as \\U = 10\\, however this can be overridden by
including a column `w_mat25` in the species parameter dataframe that
specifies the weight at which 25% of individuals are mature, which sets
\\U = \log(3) / \log(w\_{mat} / w\_{mat25}).\\

The sigmoidal function given above would strictly reach 1 only
asymptotically. For computational simplicity, any proportion smaller
than `1e-8` is set to `0`.

### Investment into reproduction

If the the energy available to a mature individual that is invested into
reproduction is not supplied via the `repro_prop` argument, it is set to
the allometric form \$\${\tt repro\\prop}(w) =
\left(\frac{w}{w\_{\tt{repro\\max}}}\right)^{m-n}.\$\$ Here \\n\\ is the
scaling exponent of the energy income rate. Hence the exponent \\m\\
determines the scaling of the investment into reproduction for mature
individuals. By default it is chosen to be \\m = 1\\ so that the rate at
which energy is invested into reproduction scales linearly with the
size. This default can be overridden by including a column `m` in the
species parameter dataframe. The sizes \\w\_{repro\\max}\\ are taken
from the `w_repro_max` column in the species parameter data frame, if it
exists, or otherwise from the `w_inf` column.

The total proportion of energy invested into reproduction of an
individual of size \\w\\ is then \$\$\psi(w) = {\tt maturity}(w){\tt
repro\\prop}(w)\$\$ In mizer edition 1, at sizes above `w_repro_max` the
value of \\\psi\\ is additionally forced to 1, so that all available
energy is invested into reproduction and growth stops. In edition 2 and
above this forcing is not applied, and \\\psi\\ is determined entirely
by the maturity ogive and the reproductive proportion.

### Reproductive efficiency

The reproductive efficiency \\\epsilon\\, i.e., the proportion of energy
allocated to reproduction that results in egg biomass, is set through
the `erepro` column in the species_params data frame. If that is not
provided, the default is set to 1 (which you will want to override). The
offspring biomass divided by the egg biomass gives the rate of egg
production, returned by
[`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md):
\$\$R\_{di} = \frac{\epsilon}{2 w\_{min}} \int N(w) E_r(w) \psi(w) \\
dw\$\$

### Density dependence

The stock-recruitment relationship is an emergent phenomenon in mizer,
with several sources of density dependence. Firstly, the amount of
energy invested into reproduction depends on the energy income of the
spawners, which is density-dependent due to competition for prey.
Secondly, the proportion of larvae that grow up to recruitment size
depends on the larval mortality, which depends on the density of
predators, and on larval growth rate, which depends on density of prey.

Finally, to encode all the density dependence in the stock-recruitment
relationship that is not already included in the other two sources of
density dependence, mizer puts the the density-independent rate of egg
production through a density-dependence function. The result is returned
by [`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md). The
name of the density-dependence function is specified by the `RDD`
argument. The default is the Beverton-Holt function
[`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md),
which requires an `R_max` column in the species_params data frame giving
the maximum egg production rate. If this column does not exist, it is
initialised to `Inf`, leading to no density-dependence. Other functions
provided by mizer are
[`RickerRDD()`](https://sizespectrum.org/mizer/reference/RickerRDD.md)
and
[`SheperdRDD()`](https://sizespectrum.org/mizer/reference/SheperdRDD.md)
and you can easily use these as models for writing your own functions.

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

## Setting resource dynamics

The `resource_dynamics` argument allows you to choose the resource
dynamics function. By default, mizer uses a semichemostat model to
describe the resource dynamics in each size class independently. This
semichemostat dynamics is implemented by the function
[`resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.md).
You can change that to use a logistic model implemented by
[`resource_logistic()`](https://sizespectrum.org/mizer/reference/resource_logistic.md)
or you can use
[`resource_constant()`](https://sizespectrum.org/mizer/reference/resource_constant.md)
which keeps the resource constant or you can write your own function.

Both the
[`resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.md)
and the
[`resource_logistic()`](https://sizespectrum.org/mizer/reference/resource_logistic.md)
dynamics are parametrised in terms of a size-dependent birth rate
\\r_R(w)\\ and a size-dependent capacity \\c_R\\. The help pages of
these functions give the details.

The `resource_rate` argument can be a vector (with the same length as
`w_full(params)`) specifying the intrinsic resource birth rate for each
size class. Alternatively it can be a single number that is used as the
coefficient in a power law: then the intrinsic birth rate \\r_R(w)\\ at
size \\w\\ is set to \$\$r_R(w) = r_R w^{n-1}.\$\$ The power-law
exponent \\n\\ is taken from the `n` argument.

The `resource_capacity` argument can be a vector specifying the
intrinsic resource carrying capacity for each size class. Alternatively
it can be a single number that is used as the coefficient in a truncated
power law: then the intrinsic carrying capacity \\c_R(w)\\ at size \\w\\
is set to \$\$c_R(w) = c_R\\ w^{-\lambda}\$\$ for all \\w\\ less than
`w_pp_cutoff` and zero for larger sizes. The power-law exponent
\\\lambda\\ is taken from the `lambda` argument.

The values for `lambda`, `n` and `w_pp_cutoff` are stored in a list in
the `resource_params` slot of the MizerParams object so that they can be
re-used automatically in the future. That list can be accessed with
[`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md).

## See also

Other functions for setting up models:
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md),
[`newSingleSpeciesParams()`](https://sizespectrum.org/mizer/reference/newSingleSpeciesParams.md),
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)

## Examples

``` r
params <- newMultispeciesParams(NS_species_params)
#> Because you have n != p, the default value for `h` is not very good.
#> Because the age at maturity is not known, I need to fall back to using
#> von Bertalanffy parameters, where available, and this is not reliable.
#> No ks column so calculating from critical feeding level.
#> Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.
#> Using f0, h, lambda, kappa and the predation kernel to calculate gamma.
```
