## Upgrading from mizer 3.3 to 3.4

The headline change is one you should not be able to feel: `MizerParams` and
`MizerSim` are now ordinary S3 objects rather than S4 objects. Every model
result is unchanged and `params@w` still works, but code that treats a model as
an S4 object needs editing, and every extension package needs a version built
for this release. Results move for a model that lets mizer fill in a `gamma` or
`f0` while an extension, external encounter or a component adds to the encounter
rate, for one whose species parameters carry an `l_mat25` that mizer had to
reject, for a calibration run that a size class holding a trace of fish was
holding up, and for a model that was rescaled with `scaleModel()` and then had
its species parameters recalculated. `summary()` of a rate array now reports
different numbers, because it used to describe sizes the species never reaches.
A model loaded with `readParams()` now has its species parameters reconciled,
which changes nothing about the model but stops values that were written
straight into the `species_params` slot from reverting later.
Everything else here is a new report or a new way of saying something, and
`info_level = 0` silences the reports.

### `MizerParams` and `MizerSim` are ordinary lists

`MizerParams` and `MizerSim` are now S3 classes. A model is a named list whose
class attribute is `"MizerParams"`, where it used to be an S4 object with formal
slots. Nothing about the model itself changed: the same slots hold the same
values under the same names, every function behaves as it did, and results are
identical to the last digit.

Slot access is unchanged too. `params@w` and `params@w <- value` still work,
through `@` methods that mizer defines for the class, so scripts written against
any earlier version keep running. What is new is that the list operators work as
well:

```r
params@w        # as before
params$w        # the same thing
params[["w"]]   # the same thing
names(params)   # the slot names
```

Models saved by an earlier version load and run as before; mizer converts them
to the new representation when it meets them, and `readParams()` and
`readSim()` return converted objects.

What no longer applies is R's S4 machinery. A model is not an S4 object, so:

| S4 code that no longer works | Use instead |
|---|---|
| `isS4(params)`, `expect_s4_class(params, "MizerParams")` | `inherits(params, "MizerParams")`, `expect_s3_class(params, "MizerParams")` |
| `slot(params, "w")`, `slot(params, "w") <- value` | `params@w`, or `params[["w"]]` |
| `slotNames("MizerParams")`, `getSlots("MizerParams")` | `names(params)` |
| `validObject(params)` | `validParams(params)` |
| `new("MizerParams")` | `emptyParams()`, or one of the `new…Params()` constructors |
| `setMethod("foo", "MizerParams", …)`, `setValidity("MizerParams", …)` | an S3 method, `foo.MizerParams` |
| `setClass("MyClass", representation(params = "MizerParams"))` | keep the model in an ordinary list element |

`setMethod()` on `"MizerParams"` is the one to watch for: it warns that there is
no definition for the class but goes on to create a method, which then never
dispatches, because the object it would dispatch on is not S4. Write
`foo.MizerParams` instead.

`is(params, "MizerParams")` and `inherits(params, "MizerParams")` both still
answer `TRUE`, as do the same tests for `"MizerSim"`, so the usual way of
checking what you were handed is unaffected.

### Extension packages need a version built for mizer 3.4

An extension package used to chain onto `MizerParams` through an S4 marker class
that mizer created at run time, and announced itself with
`registerExtension()` in its `.onLoad()`. All of that is gone: an extension now
simply prepends its own class to the S3 class vector, as in
`c("mizerShelf", "MizerParams")`. `registerExtension()` and
`registerExtensions()` have been removed, so a call to either — in a script of
your own, or in an old package's load hook — errors with `could not find
function`.

A package written against mizer 3.3 or earlier therefore has to be updated
before it can be used with mizer 3.4. If you rely on one, upgrade it and mizer
together, taking a release of the extension that names mizer 3.4; if you
maintain the package yourself, the `upgrade-extension-package` skill takes you
through the migration.

Models saved with the older package still load. The extension requirement is
recorded in the object itself, so `readParams()` and `readSim()` load the named
package and rebuild the class vector from that record.

### Saved extension objects retain their class

`saveParams()` and `saveSim()` now write the complete S3 class vector. A bare
`readRDS()` of a newly saved extension object therefore returns the extension
class as well as the base `MizerParams` or `MizerSim` class, where files written
by earlier versions deliberately contained only the base class. Code that
inspects the raw RDS file should expect the same class as the object passed to
the save helper.

Continue to use `readParams()` and `readSim()` for normal loading. Their extra
work is still needed to load required extension packages, upgrade old objects
and validate them. Validation also reconstructs the extension class of a legacy
base-class file, so existing saved models remain compatible.

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
removes it. Like `other_params()`, the new accessors show only what does not
belong to a component, so assigning a whole list through the accessor can no
longer wipe a component's entry out (#579).

Two things can break existing code:

- **Encounter contributions are now passed the current simulation time as `t`**,
  as mortality contributions already were. This makes a time-dependent extra
  encounter possible, but an existing component `encounter_fun` whose exact
  signature accepts neither `t` nor `...` will now error. Add one of those
  arguments to the function.
- **`setComponent()` refuses a component name that a free-standing contribution
  is already registered under**, where it used to take the entry over silently,
  leaving it live but invisible to the accessors. This can only affect code that
  both writes into `@other_mort` (or `@other_encounter`) by hand and then
  creates a component of the same name.

### The defaults for `gamma` and `f0` are measured on mizer's own reference state

`get_gamma_default()` works out the search volume coefficient `gamma` that gives
a species the target feeding level `f0`, by giving the species a search volume
coefficient of 1, putting a power-law resource spectrum in front of it and
measuring the energy that becomes available. `get_f0_default()` is its inverse.
Two more things used to leak into that measurement, and neither does now. (A
third, a search volume you had set by hand, was taken out of it in 3.3; see
*Defaults for `gamma` and `f0` ignore a hand-set search volume* in the 3.2 to
3.3 notes below.) Both change the species parameters of the models they
affected, and the new values are the right ones.

**The extension chain.** The measurement used `getEncounter()`, which on an
object carrying an extension class dispatches through the extension's
`projectEncounter()` method. So an extension that changes the encounter rate had
its change folded into mizer's own `gamma` — the one quantity the function goes
out of its way to derive from the species parameters alone. Both functions now
measure with `mizerEncounter()` (#577).

The consequence was not a one-off offset. `gamma` is a calculated species
parameter, so it is recomputed on every rebuild, and the recomputed `gamma` is
what the search volume is then built from — so the extension's factor was
re-applied each time: an extension that halves the search volume doubled `gamma`
on every assignment to `species_params()`. Where an extension's factor is zero
for some species — a therMizer species sitting exactly on its `temp_min`,
evaluated at the arbitrary `t = 0` the default calculation uses — the available
energy came out as zero and the calculation failed outright. An encounter
function registered with `setRateFunction()` no longer enters these defaults
either; if your model relied on that, supply `gamma` or `f0` explicitly.

**Additive encounter contributions.** The `ext_encounter` array and every
function registered with `other_encounter()`, including a component's
`encounter_fun`, are now excluded too (#586). They used to enter the call to
`mizerEncounter()` used for the measurement, combined with a search volume
coefficient of 1 for `get_gamma_default()` — where they were negligible beside
the predation encounter — but with the species' real `gamma` for
`get_f0_default()`, where they counted at full strength. The two therefore
stopped being inverses, and a model with extra food could report a substantially
higher calculated `f0` than the target its `gamma` had been derived from. A
lower recalculated `f0` is the intended reference-state value; use
`getFeedingLevel()` to inspect feeding with the extra sources present.

When no default `gamma` can be calculated the error now names the species
concerned and reports the energy measured for them. One model that used to build
now raises it: a species with no predation encounter in the reference state at
all — `interaction_resource = 0`, or a predation kernel that does not overlap
the resource — but with an additive contribution to fall back on. The `gamma` it
used to get was derived from that contribution alone and was many orders of
magnitude away from a search volume coefficient, so the model it produced was
not one to keep. Supply `gamma` explicitly for such a species.

**If you worked around any of this** by declaring the current `gamma` as a given
species parameter, you can drop the workaround and let `gamma` follow `f0`
again:

```r
given_species_params(params)$gamma <- NULL
```

The `gamma` that mizer then calculates describes the species' baseline search
volume, which is what a dynamic modulation such as a temperature scalar is meant
to modulate.

### A column dropped from an assigned species parameter table is removed

The species parameter setters read the *absence* of a column as an instruction
as well as its value. A column that the table you assign does not have is one
you no longer supply, so it is taken out of `given_species_params()`. Mizer then
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
setting that entry to `NA` (#578).

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

This is also the supported way to withdraw a species parameter column
altogether, replacing the `params@species_params[["my_col"]] <- NULL` slot
manipulation that was the only thing that worked before.

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
it, will differ (#580).

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

### The misspelled-column check runs once

The warning that a species parameter column name looks like a typo of a
standard one used to run inside every validation pass, so building a model
reported the same column many times over — 31 warnings from one
`newMultispeciesParams()` call — and every later edit of an unrelated parameter
reported it again. It now runs once, where a column enters the model, and looks
only at the columns the model does not already have (#581). The report also
goes through the mechanism that carries everything else mizer says, so
`info_level = 0` silences it; see *One report, one switch* in the 3.2 to 3.3
notes below.

Two things follow: code that built models quietly with `info_level = 0` no
longer sees the warning at all, and a test wrapping model construction in
`expect_warning()` with a count, or in `suppressWarnings()` to swallow a known
run of duplicates, is now looking at one warning rather than many. The message
itself is unchanged: `"very close to standard parameter names"`.

### `summary()` of an array covers the same sizes as `plot()`

`plot()` of a species-by-size array draws each species over its own size range,
from its `w_min` to its `w_max`. `summary()` of the same array reduced the whole
size grid, so the two described different arrays. The values outside a species'
range describe an animal that does not exist, and because a rate usually grows
with size they are also the extreme ones, so it was `Min` and `Max` that were
reported wrongly:

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
If a user reports that a summary number "changed by orders of magnitude", check
whether it came from `summary()` of a rate array before looking anywhere else —
this is much the largest such change, and the new number is the right one.
`summary(x, all.sizes = TRUE)` reproduces the old value exactly, which is the
quickest way to confirm the diagnosis.
<!-- /agent-only -->

### A size class holding no fish no longer blocks convergence

Above a size where the growth rate vanishes, a species' density decays
exponentially and `dN/dt` decays with it, so the density in the class falls
through 1e-100 and beyond while never reaching zero. `distanceSSLogN()` counted
every class with a positive density, so `log(n)` in such a class fell by the
same amount between every pair of states and its contribution to the distance
never shrank. One trace holding 3e-92 g of fish could therefore hold the
distance above any tolerance indefinitely, and the run stopped at `t_max`
reporting `converged = FALSE` (#570).

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
`matchBiomasses()`, which is a size in grams. What decides whether a state is a
fixed point is untouched: `residual_tol` is measured against a biomass drift
that integrates over every size class, cut off or not, so the cutoff cannot make
a drifting model look settled.

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

### `readParams()` reconciles the species parameters

Mizer protects the species parameters you gave it against being recalculated,
but only the ones it knows you gave, the ones in `given_species_params()`. A
value written straight into the slot with `params@species_params$h[1] <- 20` is
not among them, and the next recalculation of the species parameters -- any
change to a species parameter, and even a no-op
`species_params(params) <- species_params(params)` -- silently replaces it.
Models saved before mizer kept that record, and models built by code that
changes the slot directly, can hold many such values.

`readParams()` now calls the new `reconcileSpeciesParams()`, which records the
values a recalculation would change among the given species parameters. It
repeats that until the species parameters reproduce themselves, so a parameter
mizer derives from a hand-set one is recorded too, and the loaded model's
species parameters are a fixed point: no recalculation moves them again. It
tells you what it recorded:

```
#> The species parameter `h` holds a value that a recalculation would not
#> reproduce. I have recorded it among the given species parameters so that it
#> is not overwritten.
```

Loading a model does not change the model: no species parameter value and no
rate array is touched, only the record of where the values came from. What
changes is what happens *afterwards*. Parameters that used to revert on the
next recalculation now stay put. If you were relying on the old behaviour --
and some scripts do, without knowing it, by editing the slot and then letting a
later call undo the edit -- those parameters now keep the values the model was
loaded with.

Reaching the fixed point means recording values mizer calculated as well as
values you supplied: the `gamma` that was derived from an `h` you later
overwrote is recorded at the value the model's rate arrays were actually built
from. `calculated_species_params()` of a loaded model can therefore be smaller
than it was, and those parameters no longer respond to the parameters they were
derived from.

If you see the message on a model you did not expect to be inconsistent, the
code that built it is writing into `params@species_params` directly. Change it
to go through `species_params<-()`, or, where it has already updated the
affected rate arrays itself, through `record_given_species_params()`.

To hand one of the recorded parameters back to mizer's calculation, clear its
entry to `NA` in `given_species_params(params)`, as for any other given species
parameter.


### `scaleModel()` records the rescaled parameters as given

`scaleModel()` multiplies `R_max` by the scale factor and divides `gamma` by
it, alongside the abundances, the resource carrying capacity and the search
volume. It used to write those two values straight into the model's species
parameter table, which left them outside the record of what the user has
supplied. `given_species_params()` therefore still held the values from before
the rescaling, and the next recalculation of the species parameters -- any
change to a species parameter, and even a no-op
`species_params(params) <- species_params(params)` -- silently put them back,
undoing the rescaling with no message.

`scaleModel()` now records the rescaled `R_max` and `gamma` as given species
parameters, so they survive a later recalculation. This is the same fix
`matchGrowth()` had already received.

For most code nothing changes: a model that is rescaled and then used is
identical to before, because the values in `species_params()` were correct
immediately after the call. What changes is a model that was rescaled and then
had its species parameters recalculated. Such a model used to end up with an
`R_max` (and `gamma`) belonging to the *unscaled* model while the rest of it
was scaled, which was never what the call asked for. The new numbers are the
right ones. If that model was calibrated afterwards, the calibration was made
against the inconsistent model and is worth redoing.

`calibrateBiomass()`, `calibrateNumber()`, `matchBiomasses()` and
`matchNumbers()` all rescale through `scaleModel()` and inherit the fix.

If you have written your own version of this rescaling -- an extension package
with its own `scale...Model()` function is the likely case -- apply the same
treatment there: change the species parameters in a copy of the table and
assign it with

```r
species_params(params, recalculate = FALSE) <- sp
```

rather than writing into `params@species_params`. The `recalculate = FALSE`
records the values without rebuilding the rate arrays, which is what you want
