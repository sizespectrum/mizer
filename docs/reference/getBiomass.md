# Calculate the total biomass of each species within a size range at each time step.

Calculates the total biomass through time within user defined size
limits. By default the full size range of each species is included. With
`use_cutoff = TRUE`, the size range starts at the `biomass_cutoff`
species parameter, if it is set. You can specify minimum and maximum
weight or length range for the species. Lengths take precedence over
weights (i.e. if both min_l and min_w are supplied, only min_l will be
used).

## Usage

``` r
getBiomass(object, use_cutoff = FALSE, ..., include = "species")
```

## Arguments

- object:

  An object of class `MizerParams` or `MizerSim`.

- use_cutoff:

  If TRUE, the `biomass_cutoff` column in the species parameters is used
  as the minimum weight for each species (ignoring any size range
  arguments in `...`). If FALSE (default), the specified size range
  arguments are used, if provided, or the full size range of the species
  is used.

- ...:

  Size filters `min_w`, `max_w`, `min_l` and `max_l` as described in
  [`get_size_range_array()`](https://sizespectrum.org/mizer/reference/get_size_range_array.md).
  With a `MizerParams` object, `n`, `n_pp`, `n_other` and `t` can
  specify the state and time to measure, defaulting to the initial state
  and time zero.

- include:

  Character vector selecting `"species"` (the default), `"resource"`,
  `"components"`, a particular `"component:name"`, or an exact
  `"component:name/series"` identifier. Overlapping selections are
  included only once. Resource length limits use
  [`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md);
  resource size limits must be scalar. Components must explicitly
  declare which filters their biomass callback supports.

## Value

If called with a MizerParams object, a named vector with the biomass in
grams for each species in the model. If called with a MizerSim object,
an `ArrayTimeBySpecies` object (time x species) containing the biomass
in grams at each time step for all species. Selecting resource or
component series returns a `MizerBiomass` numeric vector for
`MizerParams`, or an `ArrayTimeBySpecies` with additional named columns
for `MizerSim`.

## Details

When `use_cutoff = TRUE`, the function checks if the `biomass_cutoff`
column exists in the species parameters. If it does, those values are
used as the minimum weight for each species. For species with NA values
in `biomass_cutoff`, the default minimum weight (smallest weight in the
model) is used.

Biomass measurement is independent of convergence assessment. The
resource is integrated over its full size grid, using the model's
quadrature scheme. Component biomass comes from the component's
registered `biomass_fun`, which receives the complete state at each
saved time. A component may return several named biomass series. Missing
measurement support or unsupported explicit filters produce NA values
with reasons in the `unavailable` attribute. Callback errors and invalid
results are reported as errors. Extended results carry `series` metadata
and a `coverage` table with `time`, `id`, `available` and `reason`
columns. Different selected series must represent different physical
sources, to avoid counting biomass twice.

## See also

Other summary functions:
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getEncounterByPrey()`](https://sizespectrum.org/mizer/reference/getEncounterByPrey.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getMortByCause()`](https://sizespectrum.org/mizer/reference/getMortByCause.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
biomass <- getBiomass(NS_sim)
biomass["1972", "Herring"]
#> [1] 218218354800
biomass <- getBiomass(NS_sim, min_w = 10, max_w = 1000)
biomass["1972", "Herring"]
#> [1] 154290431041

# If species_params contains a `biomass_cutoff`` column, it can be used
# as the minimum weight when use_cutoff = TRUE
species_params(NS_sim@params)$biomass_cutoff <- 10
biomass <- getBiomass(NS_sim, use_cutoff = TRUE)  # Uses biomass_cutoff as min_w
biomass["1972", "Herring"]
#> [1] 154290431041
```
