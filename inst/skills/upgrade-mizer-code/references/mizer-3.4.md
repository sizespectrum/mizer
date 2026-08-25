## Upgrading from mizer 3.3 to 3.4

Three things change: what `getSteadyResidual()` measures, which sizes
`summary()` of an array covers, and which size classes `distanceSSLogN()`
counts. Everything that decides whether a state is a fixed point is untouched.

### `getSteadyResidual()` measures biomass drift by default

`getSteadyResidual()` used to return a per-capita rate of change,
`(dN_i(w)/dt) / N_i(w)`, for every size class. It now returns, by default, the
contribution of each size class to the relative rate of change of its species'
*biomass*:

```r
rowSums(getSteadyResidual(params))    # (dB_i/dt) / B_i, one per species
```

The largest of those numbers, together with the resource and any other
components, is exactly what `isSteady()`, the `summary()` line of a
`MizerParams` object and `project(check_steady = TRUE)` judge a model by. The
array therefore now says *where* a model is unsteady in the same currency that
mizer uses to decide *whether* it is, and `max(abs(...))` of it is no longer a
trap.

It was a trap before. The per-capita rate of a single size class is dominated by
the fastest-relaxing cells, and those are the ones holding no mass: on
`NS_params` the largest per-capita rate is 1.2/year, in a cell holding 2e-8 of
its species' biomass, while the biomass drift is 0.014/year and the model counts
as settled. Under the second-order scheme the gap is wide enough to reverse the
ordering between a converged model and one that has just been knocked off its
steady state.

The numbers are therefore much smaller than they were, and a `plot()` of the
result has a different shape — mass-carrying sizes now stand out where fast
empty ones used to. Nothing about the model has changed. To get the old
quantity back:

```r
getSteadyResidual(params, measure = "per_capita")
```

That measure is still the right one for a question about the *structure* of a
model rather than about its drift: it shows a size class whose growth and
mortality are out of balance even where the class holds a millionth of its
species' biomass. Just do not reduce it to its maximum.

One further consequence of the new default: a size class with no fish in it is
no longer `NA`. `dN/dt` is well defined there — an empty class can be filling up
— so its contribution is reported like any other. Only a species with no biomass
at all is `NA`, for every one of its size classes. Code that passed
`na.rm = TRUE` is unaffected; code that relied on `is.na()` to find the empty
classes should test `initialN(params)` instead. With `measure = "per_capita"`
the `NA`s are exactly where they were.

The `resource` attribute of the result follows the same measure, so it too sums
to the resource's relative rate of biomass change. The `other` attribute is
unchanged: components with their own dynamics have no size grid to weight over,
so their rates stay relative ones.

<!-- agent-only -->
If a user reports that a model "became less steady" or "more steady" after
upgrading, check which number they are reading. The verdict — `isSteady()`,
`summary()`, `attr(params, "convergence")$residual` — was already
biomass-weighted in 3.3 and is unchanged. Only the diagnostic array moved. A
drop of one to three orders of magnitude in `max(abs(getSteadyResidual(params)))`
is the expected sign of the new default, not a change in the model.
<!-- /agent-only -->

### `summary()` of an array covers the same sizes as `plot()`

`plot()` of a species-by-size array draws each species over its own size range,
from its `w_min` to its `w_max`, and has done since it gained `all.sizes`.
`summary()` of the same array reduced the whole size grid, so the two described
different arrays. The values outside a species' range describe an animal that
does not exist, and because a rate usually grows with size they are also the
extreme ones, so it was `Min` and `Max` that were reported wrongly:

```r
summary(getEncounter(NS_params))$per_species[1, ]
#> before:  Species Sprat   Min 0.299   Mean 2929   Max 39573
#> now:     Species Sprat   Min 0.299   Mean 37.9   Max 240
```

40000 g/year is the encounter rate a 40 kg Sprat would have. 240 g/year is the
rate over the sizes a Sprat reaches.

Both `summary()` methods with a size dimension — `ArraySpeciesBySize` and
`ArrayTimeBySpeciesBySize` — now take `all.sizes`, defaulting to `FALSE` as
`plot()` does. Pass `all.sizes = TRUE` for the whole grid:

```r
summary(getEncounter(NS_params), all.sizes = TRUE)   # as before
```

A species with no values left in range now comes back as `NA` rather than as
the `-Inf`/`Inf` and warning that `min()` and `max()` of an empty selection
give.

`print()` and `as.data.frame()` are unchanged: they show the array as it is,
without interpreting it. Nothing about the arrays themselves changed, so any
code that indexes them directly is unaffected.

<!-- agent-only -->
If a user reports that a summary number "changed by orders of magnitude" between
3.3 and 3.4, check whether it came from `summary()` of a rate array before
looking anywhere else — this is much the largest such change, and the new number
is the right one. `summary(x, all.sizes = TRUE)` reproduces the old value
exactly, which is the quickest way to confirm the diagnosis.
<!-- /agent-only -->

### A size class holding no fish no longer blocks convergence

Above a size where the growth rate vanishes, a species' density decays
exponentially and `dN/dt` decays with it, so the density in the class falls
through 1e-100 and beyond while never reaching zero. `distanceSSLogN()` counted
every class with a positive density, so `log(n)` in such a class fell by the
same amount between every pair of states and its contribution to the distance
never shrank. One trace holding 3e-92 g of fish could therefore hold the
distance above any tolerance indefinitely, and the run stopped at `t_max`
reporting `converged = FALSE` while the biomass drift correctly reported a fixed
point.

`distanceSSLogN()` now takes a `biomass_share_cutoff`, defaulting to `1e-8`: a
size class counts only if it holds at least that share of its species' biomass.
Relevance is measured as a share of biomass rather than of density, because
density falls fifteen orders or more across a healthy spectrum for entirely good
reasons, so no density threshold could tell a dying trace from real large fish.

Nothing changes for a model without such a trace. There, every class the cutoff
removes is one that already had no density at all and was already excluded, so
the number `distanceSSLogN()` returns is unchanged to the last bit and any
`distance_tol` you have tuned keeps its meaning. What changes is the model that
has one: a run that used to reach `t_max` now converges.

Pass `biomass_share_cutoff = 0`, to `distanceSSLogN()` directly or through
`projectUntilSettled()`, for the old behaviour. Note that this is unrelated to
the `biomass_cutoff` *species parameter* used by `calibrateBiomass()` and
`matchBiomasses()`, which is a size in grams. `getSteadyResidual()` has no such
argument and needs none: weighting by biomass gives a trace the weight it
deserves without a threshold.

What decides whether a state is a fixed point is untouched: `residual_tol` is
measured against a biomass drift that integrates over every size class, cut off
or not, so the cutoff cannot make a drifting model look settled.

A run that stops at `t_max` on a state that is nonetheless a fixed point now
says so, instead of reporting only the distance:

```
#> Simulation run did not converge after 100 years. The distance function
#> returned 2.4, which is above the distance tolerance, but the state reached is
#> a fixed point: the biomasses change at only 0.0003 per year.
```

<!-- agent-only -->
If a user reports that `projectUntilSettled()` never converges on a model that
`isSteady()` calls steady, this is the first thing to check, on any mizer
version: `plot(getSteadyResidual(params, measure = "per_capita"))` and look for
a flat line at minus the mortality rate above the size where growth stops. On
3.4 and later the cutoff already handles it; before that, the remedy is to zero
the trace with `initialN(params)[initialN(params) < 1e-30] <- 0`.
<!-- /agent-only -->
