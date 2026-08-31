# Extra contributions to the mortality and encounter rates

Besides the rates that mizer calculates itself, a model can carry extra
contributions to the mortality rate and to the encounter rate. Each
contribution is an R function, registered here by name, that mizer calls
at every time step:
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md) adds
the result of every function listed in `other_mort(params)` to the
mortality rate and
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
adds the result of every function listed in `other_encounter(params)` to
the encounter rate.

## Usage

``` r
other_mort(params)

other_mort(params) <- value

other_encounter(params)

other_encounter(params) <- value
```

## Arguments

- params:

  A MizerParams object

- value:

  A named list of function names, or `NULL` to remove all contributions
  that do not belong to a component. You choose the names, but they must
  be unique and must not clash with the name of a component.

## Value

A named list with the names of the functions contributing to the rate,
excluding any that belong to a component.

## Details

Use these when the extra contribution depends on the state of the model,
as a starvation mortality or a density-dependent mortality does. A
contribution that is simply a fixed array is better set with
[`ext_mort()`](https://sizespectrum.org/mizer/reference/setExtMort.md)
or
[`ext_encounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
which is what mizer's own external mortality and external encounter use.

Assigning `NULL` to an entry removes that contribution:

    other_mort(params)[["starvation"]] <- NULL

## How your function is called

Each registered function is called as

    fun(params, n, n_pp, n_other, t, component, ...)

and has to return an array with the same dimensions as the rate it
contributes to, that is species x size. The `component` argument holds
the name under which the function is registered, so a single
implementation can serve several entries. Always give your function a
`...` argument so that it tolerates being passed arguments it does not
use.

Make sure the function depends *continuously* on the abundances, for the
reasons set out in
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md).

## Contributions belonging to a component

A component set up with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
can contribute to these rates through that function's `mort_fun` and
`encounter_fun` arguments. Such an entry belongs to its component:
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
sets it,
[`removeComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
removes it and
[`getComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
reports it. The accessors here therefore leave those entries alone.
`other_mort()` and `other_encounter()` list only the contributions that
do not belong to a component, and assigning through them preserves the
ones that do. This mirrors
[`other_params()`](https://sizespectrum.org/mizer/reference/setRateFunction.md),
which likewise hides the parameters belonging to a component, and it
makes it impossible to wipe out a component's contribution by assigning
a whole list.

## See also

[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`other_params()`](https://sizespectrum.org/mizer/reference/setRateFunction.md),
[`ext_mort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
[`ext_encounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)

## Examples

``` r
# An extra mortality that grows with the total biomass in the community.
# Like any custom rate function it has to live in the global environment
# or in a package, so that mizer can find it by name.
crowdingMort <- function(params, n, n_pp, n_other, t, component, ...) {
    biomass <- sum(n %*% (params@w * params@dw))
    # Same dimensions as the mortality rate: species x size
    0 * params@mu_b + 1e-14 * biomass
}
params <- NS_params
other_mort(params)[["crowding"]] <- "crowdingMort"
#> Error in setRateContributions(params, "other_mort", value, noun = "mortality",     arg = "mort_fun"): The entry `crowding` of `other_mort` needs to be the name of a function.
other_mort(params)
#> list()

# The contribution is included in the mortality rate
range(getMort(params) - getMort(NS_params))
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> [1] 0 0

# and can be removed again
other_mort(params)[["crowding"]] <- NULL
other_mort(params)
#> list()
```
