# Calculate the rate at which biomass of each species is fished

This yield rate is given in grams per year. It is calculated at each
time step saved in the MizerSim object.

## Usage

``` r
getYield(object)
```

## Arguments

- object:

  An object of class `MizerParams` or `MizerSim`.

## Value

If called with a `MizerParams` object, a named numeric vector with the
yield rate in grams per year for each species in the model. If called
with a `MizerSim` object, an `ArrayTimeBySpecies` object (time x
species) containing the yield rate in grams per year at each saved time
step.

## Details

The yield rate \\y_i(t)\\ for species \\i\\ at time \\t\\ is defined as
\$\$y_i(t)=\int\mu\_{f.i}(w, t)N_i(w, t)w dw\$\$ where \\\mu\_{f.i}(w,
t)\\ is the fishing mortality of an individual of species \\i\\ and
weight \\w\\ at time \\t\\ and \\N_i(w, t)\\ is the abundance density of
such individuals. The factor of \\w\\ converts the abundance density
into a biomass density and the integral aggregates the contribution from
all sizes.

The total catch in a time period from \\t_1\\ to \\t_2\\ is the integral
of the yield rate over that period: \$\$C = \int\_{t_1}^{t2}y_i(t)dt\$\$
In practice, as the yield rate is only available at the saved times, one
can only approximate this integral by averaging over the available yield
rates during the time period and multiplying by the time period. The
less the yield changes between the saved values, the more accurate this
approximation is. So the approximation can be improved by saving
simulation results at smaller intervals, using the `t_save` argument to
[`project()`](https://sizespectrum.org/mizer/reference/project.md). But
this is only a concern if abundances change quickly during the time
period of interest.

## See also

[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

Other summary functions:
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
yield <- getYield(NS_sim)
yield[c("1972", "2010"), c("Herring", "Cod")]
#> Yield rate (2 times x 2 species) [g/year] 
#>   Herring: min=3.05e+10 mean=5.52e+10 max=8e+10
#>   Cod: min=2.89e+11 mean=3.22e+11 max=3.55e+11

# Running simulation for another year, saving intermediate time steps
params <- finalParams(NS_sim)
sim <- project(params, t_save = 0.1, t_max = 1,
               t_start = 2010, progress_bar = FALSE)
# The yield rate for Herring decreases during the year
getYield(sim)[, "Herring"]
#>        2010      2010.1      2010.2      2010.3      2010.4      2010.5 
#> 30496241734 30406503779 30297907652 30173956406 30038604629 29896229522 
#>      2010.6      2010.7      2010.8      2010.9        2011 
#> 29751183405 29607379160 29468046935 29335685116 29212159486 
# We approximate the total catch in the year by averaging over the year
sum(getYield(sim)[1:10, "Herring"] / 10)
#> [1] 29947173834
```
