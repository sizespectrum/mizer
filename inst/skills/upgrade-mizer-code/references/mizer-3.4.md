## Upgrading from mizer 3.3 to 3.4

Four things change for everyone: what `getSteadyResidual()` measures, which
sizes `summary()` of an array covers, which size classes `distanceSSLogN()`
counts, and the fact that mizer now says when it invents a weight-length
relationship. For a model without components registered by `setComponent()`,
everything that decides whether a state is a fixed point is untouched; for a
model with one, the criterion no longer includes it. A fifth change concerns
anyone who assigns a species parameter table: a column missing from the table
you assign is now withdrawn from the model.
Three further sections concern extensions: repeated extension registration now
rebuilds the dynamic marker classes that reloading an extension package removed,
the defaults for `gamma` and `f0` are no longer measured through the extension
chain, and the extra contributions to the mortality and encounter rates have
gained accessors of their own.

### `getSteadyResidual()` measures biomass drift by default

`getSteadyResidual()` used to return a per-capita rate of change,
`(dN_i(w)/dt) / N_i(w)`, for every size class. It now returns, by default, the
contribution of each size class to the relative rate of change of its species'
*biomass*:

```r
rowSums(getSteadyResidual(params))    # (dB_i/dt) / B_i, one per species
```

The largest of those numbers, together with the resource, is exactly what
`isSteady()`, the `summary()` line of a
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
so their rates stay relative ones. What did change is that those rates are no
longer part of the verdict — see the next section.

<!-- agent-only -->
If a user reports that a model "became less steady" or "more steady" after
upgrading, check which number they are reading. The verdict — `isSteady()`,
`summary()`, `attr(params, "convergence")$residual` — was already
biomass-weighted in 3.3 and is unchanged. Only the diagnostic array moved. A
drop of one to three orders of magnitude in `max(abs(getSteadyResidual(params)))`
is the expected sign of the new default, not a change in the model.
<!-- /agent-only -->

### Other components are no longer part of the steady-state criterion

This section is only about models that register a component with
`setComponent()` and give it a `dynamics_fun` — an extra resource spectrum, a
detritus or carrion pool, a nutrient budget. Extensions that do this include
mizerMR, mizerReef and mizerShelf. Nothing below affects a model without one.

`isSteady()`, the `summary()` drift line, `project(check_steady = TRUE)`, the
`residual` and `attractor` entries of the `"convergence"` attribute, and the
guards in `getStability()` all used to fold in a component's rate of change
alongside the consumer and resource biomass drifts. They no longer do. They now
measure the consumers and the resource, and report the components separately.

**A model can therefore now be `isSteady()` while a component of it is
drifting.** If you have been reading `isSteady()` or `residual` as a statement
about the whole model, it never quite was one, and now it is explicitly not.
What to read instead:

```r
attr(getSteadyResidual(params), "other")   # per-component rates of change
summary(params)                            # lists them on their own lines
```

Two reasons for the change. First, a component's state is any R object at all,
so mizer knows neither what currency its entries are in nor which of its
dimensions, if any, is size. It could not form a biomass for it, and so had to
fall back on the largest per-cell relative rate — the reduction that is
dominated by the fastest-relaxing cells, which are the ones holding nothing. The
number that went into the criterion was in the wrong currency and usually far
too large.

Second, `tuneSteadyState()`, `findSteadyState(solver = "newton")` and the
stability analyses all hold the components fixed while they work. A criterion
that included them was one those functions could not satisfy: a model whose
consumers had settled to 1e-11 could report a drift of 6/year, `attractor` of
`NA`, and no way to improve either. That is what happened in issue #589.

The reports say which state variable is responsible now, rather than leaving the
number anonymous:

```
#> before: This model is not at its steady state: a biomass is changing at up to
#>         6.2 per year. [...] check `getSteadyResidual(params)` to see which
#>         species are moving.
#> now:    The consumers and the resource in this model are at their steady
#>         state, but the model as a whole is not. The component `MR` (up to 6.2
#>         per year) is also changing. Components are not included in the biomass
#>         drift above, and mizer's steady-state machinery does not settle them
```

To settle the components as well, project the model rather than tuning it:

```r
params <- findSteadyState(params, solver = "project")
# or, keeping the trajectory:
sim <- projectUntilSettled(params)
params <- finalParams(sim)
```

Both advance the components like every other state variable, and their stopping
rule still waits for them — there the components are live rather than held
fixed, so a run that stopped while one was moving would be stopping short of
something it could settle. `tuneSteadyState()` and
`findSteadyState(solver = "newton")` freeze them and cannot.

Whether a component should be able to re-enter the criterion, by declaring its
own reduction to a drift in 1/year, is issue #589.

<!-- agent-only -->
Recognise this one from the pairing: a user reports `isSteady()` disagreeing
with a residual they remember, or `attractor` of `NA` on a model whose species
residuals are at machine precision, and the model loads an extension that
registers a component (mizerMR, mizerReef, mizerShelf). In 3.3 the symptom was
the false `FALSE`; in 3.4 it is the reverse, a `TRUE` on a model with a moving
component. Both times the answer is the same: read
`attr(getSteadyResidual(params), "other")`, and settle the components with
`projectUntilSettled()` rather than `tuneSteadyState()`, which freezes them.

Do not tell the user to judge convergence from species-only residuals computed
by hand. That is the workaround that hid a 2-3/year resource drift in the report
behind #589. Mizer's own reports name the moving component; use those.
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

### Mizer says when it defaults the weight-length parameters

A species parameter data frame with no `a` or `b` column has always been given
the defaults `a = 0.01` and `b = 3`, silently. Building a model from such a
data frame now reports them among the other defaults it fills in:

```
i No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
i No `b` column so using the isometric default b = 3 in w = a l^b.
```

Nothing about the model changes; only the report is new, and it is at
`info_level = 3`, so `info_level = 0` silences it as it does the rest of the
chatter.

It is worth reading rather than silencing. `a = 0.01` with weights in grams and
lengths in cm is Fulton's condition factor K = 1, the textbook fusiform fish,
and `b = 3` is isometry. The exponent is a good default — it varies little
across species — but `a` follows body shape over roughly two orders of
magnitude, from around 0.001 for eel-like species to a few hundredths for
deep-bodied ones, so for a specific species it is a placeholder rather than an
estimate. On the twelve North Sea species shipped with mizer, whose fitted
values span `a` = 0.001 to 0.010, the default gives lengths that are 1% to 42%
short of the fitted relationship at a weight of 1 g.

That matters wherever a length reaches the model rather than a plot: the
`sigmoid_length`, `double_sigmoid_length` and `knife_edge_length` selectivity
functions convert `l50`, `l25` and `knife_edge_length` from `gear_params()`
into weights through `a` and `b`, so a defaulted relationship puts the
selectivity curve at the wrong weights and changes the dynamics. The same
applies to `min_l` and `max_l` in the summary and indicator functions, to
`getMeanLength()`, and to plots drawn with `size_axis = "l"`.

Because that case is the one that changes results, `setFishing()` reports it a
second time, where the conversion actually happens:

```
i The gear selectivity for Sprat, Cod is set by length, but `a` and `b` were
  not supplied, so the conversion to weight used mizer's defaults (a = 0.01,
  b = 3). The selectivity therefore sits at weights that are unlikely to be the
  ones you intend. Supply the weight-length parameters in the species
  parameters.
```

This one is at `info_level = 1`, so it survives the setting that silences the
routine chatter, and it is shown even when `calc_selectivity()` or
`setFishing()` is called on its own rather than through `setParams()`. It names
only the species whose selectivity is set from a length and whose `a` or `b`
mizer had to fill in, and only the parameter that was missing — supply `a` and
`b` for those species and it goes away. A gear that selects on weight
(`knife_edge`, `sigmoid_weight`) never triggers it, whatever the species
parameters say.

The test is on the selectivity function's formals: mizer hands a selectivity
function the species parameters only so that it can convert between length and
weight, so a custom `sel_func` with a `species_params` argument is covered
too.

### Repeated extension registration rebuilds missing marker classes

`registerExtension()` is designed to be called from an extension package's
`.onLoad()` hook, including when `devtools::load_all()` reloads that package in
the same session. The reload removes the S4 classes the package's namespace
held, which can take a dynamic marker class out of the registered chain.
Previously the repeated registration returned without recreating it, and the
next `coerceToExtensionClass()` call failed with base R's coercion error, `no
method or default for coercing "MizerParams" to ...`.

Repeated `registerExtension()` and `registerExtensions()` calls now rebuild the
chain's dynamic marker classes while leaving the registered chain unchanged.
The whole chain is rebuilt rather than just the class that went missing,
because R prunes a removed class from the `contains` list of its subclasses: a
marker class that used to sit outside the missing one is left parented directly
on `MizerParams`, so recreating only the missing class would leave the chain
still broken. An intact chain is inspected and left untouched, so the usual
repeated registration does no work.

The repair never installs or version-checks anything. `registerExtensions()`
loads only the namespaces of the extensions you pass it, as before, and the
repeated `registerExtension()` call touches no namespaces at all — an
extension registered earlier in the session but no longer installed cannot
make a package's `.onLoad()` fail.

### A column dropped from an assigned species parameter table is removed

The species parameter setters now read the *absence* of a column as an
instruction. A column that the table you assign does not have is one you no
longer supply, so it is taken out of `given_species_params()`. Mizer then
calculates afresh the parameters it knows how to calculate, and the ones it
does not know about leave the model:

```r
species_params(params)$gamma  <- NULL   # gamma is calculated again
species_params(params)$my_col <- NULL   # my_col is gone
```

`given_species_params(params)$… <- NULL` follows the same rule. Removing a
column is reported at `info_level` 3, with a message beginning `I have removed
the species parameter column`.

Previously neither route removed anything. Dropping a column from the table
assigned to `species_params<-()` did nothing at all: the column was restored
from the given species parameters, still recorded as given. Dropping one from
the table assigned to `given_species_params<-()` did take it out of the given
table, but left its value standing in `species_params()`, where
`calculated_species_params()` reported the user's own number as one mizer had
calculated. Nor did the removal recalculate anything, so
`given_species_params(params)$gamma <- NULL` left the previously given `gamma`
in place instead of handing it back to mizer; it is now the same instruction as
setting that entry to `NA`.

Two things to watch for in existing code:

- **Assigning a table with only some of the model's columns now withdraws the
  rest.** This used to be a way of updating a few columns and leaving the others
  alone. Edit the table you get from `species_params(params)` instead of
  building a new one. The scope for surprise is limited, because
  `species_params<-()` has always validated what it is given: a table without
  `species` and one of `w_inf`, `w_max` or `w_repro_max` is an error, not a
  partial update.
- **A rebuild now happens where it did not before.** Removing a column from the
  given species parameters goes through `setParams()`, so a parameter that was
  derived from the withdrawn one moves to its calculated value. That is the
  point of the change, but it means numbers can shift where previously nothing
  did.

For extension packages this is the supported way to withdraw a species
parameter column when the user switches the extension off, replacing the
`params@species_params[["my_col"]] <- NULL` slot manipulation that was the only
thing that worked before.

### Defaults for `gamma` and `f0` ignore the extension chain

`get_gamma_default()` works out the search volume coefficient `gamma` that gives
a species the target feeding level `f0`. It does that by giving the species a
search volume coefficient of 1, putting a power-law spectrum in front of it and
measuring the energy that becomes available. It used to measure with
`getEncounter()`, which on an object carrying an extension marker class
dispatches through the extension's `projectEncounter()` method. So an extension
that changes the encounter rate had its change folded into mizer's own `gamma` —
the one quantity the function goes out of its way to derive from the species
parameters alone. `get_f0_default()`, the inverse, had the same problem. Both
now measure with `mizerEncounter()`.

The consequence was not a one-off offset. `gamma` is a calculated species
parameter, so it is recomputed on every rebuild, and the recomputed `gamma` is
what the search volume is then built from — so the extension's factor was
re-applied each time:

```r
species_params(params)$gamma
#> Doubled on every assignment to species_params(), on an extension that
#> halves the search volume; now stays put.
```

Where an extension's factor is zero for some species — a therMizer species
sitting exactly on its `temp_min`, evaluated at the arbitrary `t = 0` the
default calculation uses — the available energy came out as zero and the
calculation failed outright with `Could not calculate gamma.`

The same now applies to an encounter function registered with
`setRateFunction()`: it no longer enters these two defaults. If your model
relied on that, supply `gamma` (or `f0`) explicitly in the species parameters
instead.

**If you are an extension author** who worked around this by declaring the
current `gamma` as a given species parameter, you can drop the workaround and
let `gamma` follow `f0` again:

```r
given_species_params(params)$gamma <- NULL
```

The `gamma` that mizer then calculates describes the species' baseline search
volume, which is what a dynamic modulation such as a temperature scalar is
meant to modulate. Models that carried a `gamma` inflated (or deflated) by the
extension's factor will get different numbers from this release, and those
numbers are the right ones.

### `other_mort()` and `other_encounter()` register contributions that have no component

`getMort()` and `getEncounter()` add the result of every function listed in
`params@other_mort` and `params@other_encounter`. Until now the only exported
way to write into those lists was `setComponent()`, which needs a
`dynamics_fun` and an `initial_value`, so an extension adding a term that
depends on the model state but keeps no state of its own — a starvation or
senescence mortality — had to assign into the slot directly. That still works,
but the supported way is now

```r
other_mort(params)[["starvation"]] <- "starvMort"
other_encounter(params)[["scavenging"]] <- "scavengingEncounter"
```

which checks that the name really is a function and that it does not collide
with a component's. Names must be unique, and assigning `NULL` to an entry
removes it.

Encounter contributions are now passed the current simulation time as `t`, as
mortality contributions already were. This makes a time-dependent extra
encounter possible. An existing component `encounter_fun` whose exact signature
accepts neither `t` nor `...` will now error; add one of those arguments to the
function.

Like `other_params()`, the new accessors show only what does not belong to a
component. An entry created by `setComponent(mort_fun = )` stays the property of
its component: `getComponent()` reports it, `removeComponent()` removes it, and
`other_mort()` does not list it — so assigning a whole list through the
accessor can no longer wipe it out.

One case that used to pass now stops. `setComponent()` refuses a component name
that a free-standing contribution is already registered under:

```
#> Error: There is already a rate contribution registered under the name
#> 'starvation'. Remove it with `other_mort()` or `other_encounter()` first, or
#> choose a different component name.
```

Previously the component silently took the entry over, which left it live but
invisible to the accessors. This can only affect code that both writes into
`@other_mort` (or `@other_encounter`) by hand and then creates a component of
the same name.

### A misspelled species parameter column is reported once, and quietly

Mizer flags a column name that looks like a typo of a standard species
parameter. It used to run that check inside every validation pass, so building
a model reported the same column many times over — 31 warnings from one
`newMultispeciesParams()` call — and every later edit of an unrelated parameter
reported it again. The check now runs once, where a column enters the model, and
looks only at the columns the model does not already have.

Two things to watch for in existing code:

- **The report follows `info_level`.** It is now raised through the same
  mechanism as everything else mizer says while it sets up a model, so
  `info_level = 0` and `options(mizer_info_level = 0)` silence it. Code that
  built models quietly and relied on a typo still being flagged no longer gets
  the warning. Build at the default `info_level` when you want the check.
- **A test that counted warnings will see a different number.** Anything
  wrapping model construction in `expect_warning()` with a count, or in
  `suppressWarnings()` to swallow a known run of duplicates, is now looking at
  one warning rather than many.

The message itself is unchanged: `"very close to standard parameter names"`.

### An invalid `w_mat25` is replaced by its default

`validSpeciesParams()` rejects a `w_mat25` that is not smaller than `w_mat`. It
used to say that it had corrected that *by setting it to NA*, which described a
step the user never saw: `setReproduction()` fills the default in straight
afterwards. It now says it is
`"marking it as missing so that its default will"` be used.

Where the species parameters also carry `l_mat25`, this is a change in results
and not only in wording. The length used to survive the correction and the
length-to-weight conversion put the rejected value straight back, a rounding
error below `w_mat` — small enough to pass the `w_mat25 < w_mat` assertion in
`setReproduction()`, and large enough to collapse the maturity ogive to a step
function. Such a model now gets the intended default,
`w_mat / 3^(1/10)`, and a smooth ogive. Maturity, and everything downstream of
it, will differ.
