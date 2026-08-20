## Upgrading from mizer 3.2 to 3.3

Most of the changes in this release are corrections. Results move only for
models that had opted in to second-order bin-averaging or to the `van_leer`
flux, that set `min_w` below the default, that specify sizes as lengths, that
change the resource power law after constructing the model, or that were brought
to steady state with `steadyNewton()` while their consumers were satiated. The
one change to an interface is in the spectrum plots.

### `biomass` and `per_log_size` replace `power`

`plotSpectra()`, `plotSpectra2()`, `plotCDF()`, `plotCDF2()` and `animate()`
now describe the plotted quantity with two independent arguments: `biomass`
chooses a biomass density rather than a number density, and the new
`per_log_size` chooses a density with respect to logarithmic size rather than
with respect to size. The `power` of the weight multiplying the number density
is the sum of the two:

| | `per_log_size = FALSE` | `per_log_size = TRUE` |
|---|---|---|
| `biomass = FALSE` | `power = 0` | `power = 1` |
| `biomass = TRUE`  | `power = 1` | `power = 2` |

`power` still works and is still the only way to ask for a power that is not
the sum of the two flags, so calls that pass only `power` are unaffected. Two
things change:

- Passing `power` together with `biomass` used to ignore `biomass` silently.
  Now the two must agree, or you get an error (#501). Where they do agree the
  call is honoured: `plotSpectra(sim, power = 1, biomass = FALSE)` used to plot
  the biomass density, and now plots the number density with respect to
  logarithmic size — the same numbers, but labelled correctly and, with
  `size_axis = "l"`, converted to a length axis with the logarithmic Jacobian
  rather than the density one. If you meant the biomass density, drop the
  `biomass` argument. The same applies to `plotlySpectra()`, `plotlyCDF()`,
  `plotlySpectra2()` and `plotlyCDF2()`, which passed `power` on internally and
  so ignored `biomass` even when you gave only `biomass`: those calls now plot
  what they were asked for.
- `plotCDF()` and `plotCDF2()` do not accept `per_log_size`, because
  integrating a density over size gives the same cumulative quantity either
  way. Use `biomass` on its own there.

### Arrays say what kind of value they hold

Mizer arrays now carry a `type` attribute saying what kind of quantity their
values are: `"value"` (the default) for a rate or an amount, `"density"` for an
amount per gram of body weight, `"proportion"` for a fraction. Two kinds of
plotting behaviour follow from it, and both used to be decided some other way.

**Densities.** Plotting a density against a length axis (`size_axis = "l"`) has
to multiply the values by a Jacobian, because a density per gram is not a
density per centimetre. mizer used to decide which arrays those were by looking
at their metadata strings, treating an array as a density if it was named
`"Number density"` or had units `"1/g"`. For mizer's own number spectra —
`initialN()`, `N()`, `finalN()`, `NResource()`, `resource_capacity()` — nothing
changes; they were recognised before and are tagged now. What changes is
`getFluxGradient()`: it is a rate of change of a number density, with units
`g^-1/year`, and neither of the old string tests recognised it, so on a length
axis its values were left as densities per gram and were mislabelled as such.
They are now converted with the `dw/dl = b w / l` Jacobian and labelled
`cm^-1/year`. The new curve is the right one; if you were reading values off the
old one, they were per gram plotted against length.

**Proportions.** `getFeedingLevel()`, `getCriticalFeedingLevel()`, `maturity()`,
`repro_prop()`, `psi()` and `resource_level()` now declare themselves
proportions, and a plot of one shows the whole of the interval from 0 to 1 on a
linear y axis, so the value can be read against the scale it belongs to. Three
consequences:

- `plot(getFeedingLevel(params))` and the other array plots gain that y range,
  where they used to fit the axis to the data. This is the range
  `plotFeedingLevel()` has always shown, so the dedicated function and the array
  plot now agree.
- `plot(resource_level(params))` gets a linear y axis instead of a logarithmic
  one. Pass `log_y = TRUE` to get the old axis back; any explicit `log_y` or
  `log` you already pass is respected.
- The range is only ever *widened* to include the data, never narrowed to the
  interval from 0 to 1. So `plotFeedingLevel(include_critical = TRUE)` now shows
  a critical feeding level above 1, which the old fixed window drew off the top
  of the plot. Nothing is ever hidden, and an explicit `ylim` still wins.

**Declaring it yourself.** An array of your own is taken to be a density or a
proportion only if you say so, by passing `type` to the array constructor. If
you do not pass it, the old string tests still run as a fallback, so existing
code that named an array `"Number density"` or gave it units `"1/g"` keeps
working, and arrays saved by earlier versions keep working when they are loaded.

Extension packages that called the unexported plotting helpers directly should
note that `plotComparisonDataFrame()` and the internal `animate_plotly()` take a
single `density_wrt` argument in place of `spectrum_power` and
`spectrum_per_log_size`, and that the internal `array_spectrum_power()` is gone.
The `power`-based interface of `plotSpectra()` and friends is unchanged.

### The total is summed on the axis it is plotted against

A total can only be formed once every line sits on the same coordinate. On a
weight axis they do: every species shares the model's weight grid. On a length
axis they do not, because each species — and now the resource — converts weight
to length with its own allometric relationship, so at a given length the lines
sit at different weights. That is why the `Total` line used to be dropped from
length-based plots.

It is now summed *after* the conversion, at equal length rather than at equal
weight, interpolating each line onto the union of all the size coordinates
(logarithmically in size, with a line contributing nothing outside its own
range). Where the lines already share a grid — always on a weight axis, and on
a length axis whenever the weight-length parameters agree — the union is that
grid and the interpolation reproduces the values exactly. **The weight-axis
total is unchanged**, for every power.

`plotSpectra2()` and `plotSpectraRelative()` are fixed by the same change.
They used to convert the size axis after assembling the two spectra, so the
total they had been handed — already summed at equal weight — reached the
conversion with no species to convert it by and was silently dropped. They now
let `plotSpectra()` do the conversion, so the total they receive is the total
on the axis being plotted.

That also settles what `ylim` does there. `plotSpectra()` applies it both as
the axis limits and as a filter on the data, with a hard floor at 1e-20.
`plotSpectra2()` could not do the same on a length axis, because the values it
was filtering were a Jacobian away from the ones the limits described, so it
skipped the filter. Now that it converts first, the filter applies as it does
everywhere else. The plot is unchanged — the axis limits hid those points
anyway — but `return_data = TRUE` no longer hands back values outside the
limits.

One thing does change on the weight axis. `total = TRUE` now means the same
thing everywhere: **the total of everything the object holds**, whatever is
drawn. `plotSpectra()` always worked that way — the resource and every species
count, whether or not the resource is shown and whichever species were
selected — and it still does. The array plots did not: `plot(<array>,
total = TRUE)` summed only the species selected for display, and only the sizes
inside each species' own range. It now sums the whole array, so the total no
longer moves when you change `species`, `all.sizes` or `background`, and a plot
of two species can be read against the community total. If you were relying on
the total of a selection, sum the selection yourself.

### Length and weight parameters follow the one you gave last

A size can be given either as a weight (`w_mat`, `w_max`, …) or as the length it
converts to (`l_mat`, `l_max`, …). Mizer used to derive the weight from the
length whenever both were present, so on a model specified by lengths a weight
could not be set at all: the value you assigned was replaced on the spot by the
one calculated from the unchanged length.

Both now follow one rule: **the one you gave last wins, and if you gave both at
the same time the weight wins.** The other is set to match, so the two never
disagree, and mizer warns, naming the species, when it changes a length to match
a weight it disagrees with.

```r
params <- newMultispeciesParams(sp)   # sp specifies l_mat, a and b

# Used to be silently undone, now it takes effect and l_mat follows
species_params(params)$w_mat[1] <- 100
```

The rule is applied when a data frame is assigned into a model, which is when
mizer can tell which values changed. A data frame you have taken out of a model
and are editing on its own is left exactly as you write it — the conversions,
checks and warnings happen on assignment. One that was never in a model, for
example one passed to `validSpeciesParams()`, carries no such history, so a
length and a weight that disagree there count as given at the same time and the
weight wins.

`species_params<-()` and `given_species_params<-()` apply the rule identically
(#490). Previously only `species_params<-()` did, so the same edit made through
`given_species_params<-()` was discarded — a maturity weight differing by up to
73% — and the given species parameters were left permanently inconsistent, which
made mizer repeat the "not consistent" warning at every later parameter change.
If you worked around this by setting the weight and the length together, you can
now set either one on its own.

### Species parameter setters distinguish edits from declarations

The setters now express two different intentions:

| Setter | Meaning |
|---|---|
| `species_params<-()` | Edit the complete table. Mizer detects and records only entries whose values changed. |
| `given_species_params<-()` | Declare the authoritative user input. Every non-`NA` entry is given, even when equal to the current calculated value; `NA` or a removed column hands it back to mizer. |

This makes it possible to protect a calculated value without changing the
current model:

```r
given_species_params(params)$q <- species_params(params)$q
```

The setters rebuild through `setParams()` only when the model can change.
Provenance-only changes, observations, direct-runtime parameters and unrelated
custom columns on a base `MizerParams` object are stored without rebuilding all
rate arrays. Dependent parameters, demotions to calculated values, arguments of
the active predation kernel and unknown columns on extension objects retain the
conservative rebuild path. Call `setParams()` explicitly if the intention is to
repair an object after direct slot manipulation.

`given_species_params<-()` also reports instructions that cannot take effect;
`species_params<-()` stays quiet. It warns when a parameter is overruled by
another given parameter, feeds a rate array set by hand, or belongs in
`gear_params()`. Clearing an actually given value to `NA` counts as a change;
adding an all-`NA` column does not. A frozen resource array is handled the same
way. These warnings expose an existing no-op rather than changing the model:
reset the named array to return it to parameter control, set the array directly,
or use `options(mizer_info_level = 0)` when the warning is not wanted (#489).

Finally, defaults added while validating a `species_params<-()` assignment are
no longer mistaken for values the user supplied. In particular, an unrelated
edit no longer freezes newly filled `a` and `b` values. Such values can now move
again when their inputs change and appear in `calculated_species_params()`
rather than `given_species_params()`. Set a value explicitly if it should remain
fixed (#496).

### One report, one switch

Nearly everything mizer says while building or changing a model now goes
through the same mechanism, including the reports in `steady()`,
`projectToSteady()`, `validParams()`, `setInteraction()`,
`setReproduction()`, `setResource()`, `newTraitParams()`,
`newSingleSpeciesParams()` and `plotYieldObservedVsModel()`. Two consequences
for existing code:

- **`info_level = 0` now means silence.** Reports that were plain `message()`
  calls ignored `info_level` altogether and appeared anyway; they no longer do.
  If your code relied on seeing one of them, drop the `info_level = 0`.
- **Reports are collected and given at the end of the call**, one message and
  one warning rather than a stream. A test doing `expect_message()` on an
  individual report inside a longer call may need adjusting, and the text now
  arrives with any others in the same message.

### `$` on a parameter table no longer partially matches

`$` on a `species_params` or `gear_params` table now matches column names
exactly. Partial matching was silently returning the wrong parameter: in a model
without the length-weight parameters `a` and `b`,

```r
species_params(NS_params)$a   # used to return the `alpha` column
species_params(NS_params)$b   # used to return the `beta` column
```

complete with per-species names, so code converting weights to lengths got the
assimilation efficiency and the preferred predator/prey mass ratio instead.
Writing was never partially matched (`sp$b <- 3` always created a new column
`b`), so reads and writes disagreed about what `$b` meant.

A name that is not a column now gives `NULL`. If the name would have partially
matched a single column, you also get a warning naming that column. So
`is.null(species_params(params)$foo)` is now a reliable test for whether a
parameter is present, and any code that was relying on the abbreviation should
spell the column out in full (#487).

### Quadrature fixes under `second_order_w`

Two diagnostics were applying the prey-bin quadrature twice when second-order
bin-averaging was switched on with `second_order_w()`, and are now consistent
with `getEncounter()`:

- **`getDiet(proportion = FALSE)`** was uniformly too large by a factor
  `(1 + beta) / 2`, where `beta` is the grid ratio — 9.7% for `NS_params`.
  Summing the diet over prey now reproduces
  `getEncounter() * (1 - getFeedingLevel())` under both schemes (#474).
  `getDiet(proportion = TRUE)`, the default, was unaffected: the factor was
  uniform and divided out.
- **`getTrophicLevel()`** built its numerator and denominator from different
  quadratures, so reported trophic levels were off by up to 0.06. A predator
  whose prey all have trophic level 1 now comes out at exactly 2 under both
  schemes (#474).

Models on the default (first-order) scheme are unchanged. If you have published
absolute diet values or trophic levels computed under `second_order_w`, they need
recomputing.

### `getProportionOfLargeFish()` on a `MizerParams` object was wrong

The `MizerParams` method multiplied the species x size abundance array by the
vector of weights, which R recycles down the columns of the array rather than
along the size axis, so every species but the first was weighted by the wrong
sizes. Only the `MizerParams` method was affected; the `MizerSim` method was
always right, and the two now agree when applied to the same state (#494). Any
Large Fish Index computed from a `MizerParams` object in a model with more than
one species needs recomputing.

### `getN()` respects the quadrature scheme at the ends of the size range

`getN(params, min_w = ...)` now bin-averages the size-range window when
second-order bin-averaging is switched on with `second_order_w()`, so the bin
straddling `min_w` or `max_w` contributes only partially — as `getBiomass()`
already did. Numbers over a restricted size range therefore change slightly
under `second_order_w`; over the full size range, and on the default
first-order scheme, nothing changes (#494).

### The calibration and matching functions respect the quadrature scheme

`calibrateBiomass()`, `calibrateNumber()`, `matchNumbers()`,
`plotBiomassObservedVsModel()` and `plotYieldObservedVsModel()` each wrote out
their own sum over the size grid rather than using mizer's size integral, so
they stayed on the first-order quadrature and cut the size range at a bin
boundary even in a model with second-order bin-averaging switched on with
`second_order_w()`. In such a model the calibration functions left the
abundances at values that disagreed with the `getBiomass()` or `getN()` you
would check them against, and a species matched to its observed biomass was
then plotted off the 1:1 line. All five now integrate the same way
`getBiomass()`, `getN()` and `getYield()` do, so a matched species really does
come out at its observation. On the default first-order scheme nothing changes
(#504, #529).

`plotYieldObservedVsModel()` is the one where the size of the error matters.
Its model yields were 10-20% below `getYield()` for a model on the
second-order scheme, and the total relative error in the plot caption is
computed from them, so it told you the model under-predicted the yields when it
did not. If you have read a yield calibration off that plot under
`second_order_w()`, re-read it.

`matchNumbers()` also gains the empty-selection guard that `matchBiomasses()`
already had. Its own guard could never fire, so when it had nothing to match —
no `number_observed` values, or none for the species you asked for — it left the
abundances alone but still called `setBevertonHolt()`, updated `time_modified`
and announced that it had moved the model off its steady state. It now returns
the model unchanged, as `matchBiomasses()` always did. Code that relied on the
incidental re-tuning of the reproduction parameters should call
`setBevertonHolt()` itself.

### `w_min` survives a rebuild of the species parameters

`w_min` is now part of `given_species_params`, so the `min_w` argument to
`newMultispeciesParams()` and `emptyParams()` is preserved across any operation
that rebuilds the species parameters. Previously a `given_species_params<-`
round-trip silently reset `w_min` to 0.001 when `min_w` was smaller, and emitted
a spurious warning when it was larger (#460). Code that set a small `min_w` and
worked around the reset — or that unknowingly ran on the reset grid — now gets
the size grid it asked for, and results change accordingly.

### The `match…()` functions announce that they broke the steady state

`matchBiomasses()`, `matchNumbers()` and `matchGrowth()` now
report that they have moved the model off its steady state. This is a message,
so `info_level = 0` or `options(mizer_info_level = 0)` silences it, as does the
`info_level` argument, which `matchGrowth()` gains and which `matchBiomasses()`
and `matchNumbers()` previously accepted but ignored.

The `calibrate…()` functions and `scaleModel()` say nothing, because they do not
break the steady state: they apply one overall scaling factor, which is an exact
symmetry of the model. If your workflow re-ran `steady()` after every
`calibrate…()` step, that step was never necessary.

### `summary()` reports the steady state

`summary()` of a `MizerParams` object has a new block:

```
Steady state:
	biomass drift:	3.2e-05 /year	(at steady state)
```

Code that parses the output of `summary()` by line position needs updating. The
same number is available directly as `getSteadyResidual()`.

### The convergence attribute gained a `residual` field

The `"convergence"` attribute attached by `projectToSteady()` and `steady()` has
a new `residual` entry, so `expect_named()` or `names()` checks on it need
updating. It reports how far the state reached actually is from a fixed point,
which the existing `distance` field only approximates — `distance` compares two
states `t_per` apart on whatever scale the distance function uses.

When the two disagree, `steady()` now appends to its convergence message:

```
#> Convergence was achieved in 12 years. A biomass is still changing at up to
#> 0.42 per year. Reduce `tol` to converge further.
```

This is the case where the relative-RDI criterion is satisfied while the spectra
are still moving. It is a message, not a warning, because convergence at the
`tol` you asked for did happen.

### `getStability()` checks that it was given a steady state

Both `getStability()` and `getLimitCycleSim()` linearise the dynamics *at*
`initialN(params)`. If that state is not a fixed point, the eigenvalues describe
the neighbourhood of a point the model is not sitting at and the verdict on
stability is meaningless. Both now warn in that case. Run `steadyNewton()` first,
or silence with `options(mizer_info_level = 0)` if you know what you are doing.

### `steady()` converges under the `van_leer` flux scheme

On a model whose `second_order_w()` selects the `"van_leer"` flux, `steady()`
used to fall into a limit cycle instead of converging: the flux limiter weights
flipped from one cell to the next between iterations, and the iteration chased
itself. The limiter is now relaxed with an exponential moving average, and the
run converges (#522).

Code that worked around this — a `steady()` call wrapped in `try()`, a hand-set
`t_max`, a fall-back to the default upwind flux, or a `steadyNewton()`
substituted for `steady()` — is no longer needed. The steady state it now
reaches is the one the `van_leer` discretisation actually has, so it differs
from the upwind steady state the workaround was settling on; recalibrate rather
than treat the difference as a regression.

### `steadyNewton()` solves for the resource

`steadyNewton()`'s analytic substitution for the semichemostat resource assumed
that consumer feeding levels were fixed while the resource adjusted, which is not
self-consistent once consumers are satiated: the resource density and the feeding
level it produces determine each other. The resource is now carried among the
solver's unknowns, so the two are updated together (#521).

The fixed point this converges on is the correct one, so **steady states found
with `steadyNewton()` on a model with satiated consumers move**, and anything
downstream of them — `getStability()`'s spectral radius, `getLimitCycleSim()`'s
period, a `plotBifurcation()` diagram — moves with them. Models whose consumers
are far from satiation are unaffected. `getStability()`'s quasi-static
approximation gained a fixed iteration for the same reason, which also makes its
numerical Jacobian smoother; small changes in the reported eigenvalues are
expected.

### `compareParams()` compares small parameters properly

`compareParams()` now uses a relative tolerance for species parameters, so
small-magnitude parameters such as `gamma` (~1e-8) are no longer treated as equal
when they differ by up to ~10%. Comparisons that previously reported two models
as identical may now report differences — those differences were always there.

### Resource scalars refresh calculated `gamma` and `q`

The resource power law is also the reference spectrum used to calculate search
volume parameters. Changing `lambda` through `resource_params<-()` or
`setResource()` now recalculates every `q` and `gamma` entry that mizer owns;
changing `kappa` recalculates every mizer-owned `gamma`. A value you supplied
explicitly remains protected, including when only some species in a column were
given (#497).

Previously the resource capacity was rebuilt but the calculated species
parameters and `search_vol` were left at the values for the old resource:

```r
params <- newMultispeciesParams(sp)
resource_params(params)$lambda <- 2.2

# These now follow the new lambda automatically
species_params(params)$q
species_params(params)$gamma
```

If existing code deliberately wanted to keep the old values, record them as
given before changing the resource:

```r
given <- given_species_params(params)
given$q <- species_params(params)$q
given$gamma <- species_params(params)$gamma
given_species_params(params) <- given
resource_params(params)$lambda <- 2.2
```

### Defaults for `gamma` and `f0` ignore a hand-set search volume

`get_gamma_default()` works out how much energy is available to a predator by
giving it a search volume coefficient of 1. It used to obtain that search volume
by calling `setSearchVolume()`, which refuses to recalculate a `search_vol`
array you have set by hand — so mizer's own internal call was blocked along with
yours, and the available energy was measured with *your* array. The resulting
`gamma` was wrong by whatever factor separated your array from the unit-gamma
one, which in a realistic model is many orders of magnitude. `get_f0_default()`,
the inverse, had the same problem. Both now build the search volume they need
directly from the species parameters (#488).

```r
sv <- search_vol(params)
search_vol(params) <- sv * 10          # freeze the search volume

given_species_params(params)$gamma <- NA   # ask mizer to recalculate gamma

species_params(params)$gamma
#> Used to come back ~1e9 times too large; now the same value you would
#> get without the frozen search volume.
```

If you have a model in which you set `search_vol` by hand and then let mizer
fill in a missing `gamma` or `f0`, that model's species parameters were wrong
and change with this release. Note that the recalculated `gamma` still has no
effect on the model while the search volume stays frozen — mizer now warns you
about that separately, see "Species parameter setters distinguish edits from
declarations" above.
Call `setSearchVolume(params, reset = TRUE)` to put the search volume back under
the control of the species parameters.

### `f0` is always validated

Every non-missing target feeding level `f0` must now be finite and in the
interval `[0, 1)`, whether or not a search-volume coefficient `gamma` is also
supplied. Previously `f0 = 1` divided by zero when mizer calculated `gamma`,
silently creating an infinite `gamma` and a non-finite `search_vol`; values
above 1 created negative search volumes. If `gamma` was supplied explicitly,
the same invalid `f0` could instead be accepted and ignored.

If code now errors here, choose a physically attainable feeding level below 1.
When `gamma` is the parameter you intend to control, omit the `f0` value or use
a valid value; `gamma` will still take precedence (#517).

### `setExtMort()` warns when `z0pre` or `z0exp` is ignored

The `z0pre` and `z0exp` arguments of `setExtMort()` are used only to calculate
values of the `z0` species parameter that are not present in
`given_species_params()`. If `z0` is given for every species, calls such as

```r
params <- setExtMort(params, z0pre = 2)
params <- setParams(params, z0exp = -0.25)
```

were accepted but changed nothing. They now warn that the arguments were
ignored. `reset = TRUE` does not make them applicable: it hands the
external-mortality array back to the species parameters but does not remove the
given `z0` values.

Set `z0` explicitly when changing an existing model:

```r
given_species_params(params)$z0 <- 2 * species_params(params)$w_inf^(-0.25)
```

The arguments still work wherever `z0` is not given. A `z0` value present only
in `species_params()` is the cached result of an earlier calculation and is
recalculated. If either argument was supplied explicitly, the newly calculated
`z0` values are recorded in `given_species_params()` so that they survive later
rebuilds. Values calculated from the default `z0pre = 0.6` and
`z0exp = n - 1` remain calculated parameters and are not recorded there (#493).

### `setParams()` rejects arguments it does not use

`setParams()` passes its `...` on to the rate setters, each of which declares
its own `...` as unused. Any argument that none of them recognises was
therefore accepted and ignored without a word. It is now an error, and the
error says where the argument belongs when it belongs somewhere:

```r
setParams(params, metabolic = 99)        # was: silently ignored
setParams(params, resource_rate = 5)     # was: silently ignored
```

The resource case is the one most likely to have bitten: `setParams()` never
called `setResource()`, so no resource argument ever reached the model, and the
deprecation warnings for `setResource(r_pp)` and `setResource(kappa)` used to
recommend `setParams(resource_rate)` and `setParams(resource_capacity)`, which
do nothing. Use `setResource()` for all of these:

```r
params <- setResource(params, resource_rate = 5)
```

Likewise `gear_params` goes to `gear_params<-()`, and `second_order_w` and
`use_predation_diffusion` to their own assignment functions. `setResource()`
itself now applies the same check to its own `...`, so a misspelled resource
argument errors there too.

If your code errors here, the argument was having no effect before, so removing
it changes nothing; moving it to the right function changes the model, and that
is the change you had intended all along.

Two related tidy-ups: `reset` is now a documented argument of `setParams()`
(it was already forwarded through `...`, undocumented) and still thaws every
rate array that `setParams()` sets; and `setExtDiffusion()` is now listed
among the setters that `setParams()` calls, which it always did.

<!-- agent-only -->

Diagnostic: if a user reports that a resource change "did not take", check
whether they went through `setParams()`. Before this release the call was
accepted, so there is no error in their logs and the model simply kept its old
resource. `resource_rate(params)` before and after their call is the quickest
confirmation.

<!-- /agent-only -->

### One name for each stored rate array

Seventeen accessors that read a parameter or rate array back out of a `MizerParams` object had
two interchangeable names. The bare name is now the one to use — it is the one
that also has a replacement function, so the pair reads the same way in both
directions (`catchability(params)` and `catchability(params) <- value`,
`reproduction_level(params)` and `reproduction_level(params) <- value`). The
`get`-prefixed names are superseded:

| Superseded | Use instead |
|---|---|
| `getCatchability()` | `catchability()` |
| `getSelectivity()` | `selectivity()` |
| `getInitialEffort()` | `initial_effort()` |
| `getInteraction()` | `interaction_matrix()` |
| `getResourceDynamics()` | `resource_dynamics()` |
| `getResourceLevel()` | `resource_level()` |
| `getResourceRate()` | `resource_rate()` |
| `getResourceCapacity()` | `resource_capacity()` |
| `getPredKernel()` | `pred_kernel()` |
| `getSearchVolume()` | `search_vol()` |
| `getMaxIntakeRate()` | `intake_max()` |
| `getMetabolicRate()` | `metab()` |
| `getExtMort()` | `ext_mort()` |
| `getExtEncounter()` | `ext_encounter()` |
| `getMaturityProportion()` | `maturity()` |
| `getReproductionProportion()` | `repro_prop()` |
| `getReproductionLevel()` | `reproduction_level()` |

Nothing breaks: the old names are kept as plain aliases of the new ones. They
do not warn, they will not be removed, and they return exactly the same value,
so old code and old scripts keep running untouched. Renaming is a search and
replace whenever you next touch the code.

The `get` prefix now means one thing — a function that *calculates* something
from the current state of a model, like `getEncounter()`, `getFMort()` or
`getBiomass()`. The functions above only hand back a value that is already
stored in the object.

<!-- agent-only -->

The `get` forms are no longer S3 generics of their own; each is now bound to the
same function object as the bare name (`getExtMort <- ext_mort`). Dispatch still
works for a custom class, but on the bare name, so an extension that defined
`getExtMort.MyClass` or `getInteraction.MizerParams` must rename its method to
`ext_mort.MyClass` or `interaction_matrix.MizerParams`. A method on the bare
name has always worked and keeps working through both names.

Because the aliases no longer warn, a user running old code sees no signal at
all. Do not tell them their code is about to break; it is not. Suggest the new
name when you are already editing the line.

When a user's `bin_average` diagnostic disagrees with the rate functions, check
whether they reached for `pred_kernel()` (formerly `getPredKernel()`) rather
than `encounter_kernel()` — the rename does not change that distinction, but it
makes the two names look more alike than they used to.

<!-- /agent-only -->

### `yield_observed` belongs to the gear parameters

`plotYieldObservedVsModel()` now takes the observed yield from the
`yield_observed` column of `gear_params()`, where the yield is given for each
gear-species pair and the plot adds it up over the gears:

```r
gear_params(params)["Cod, Otter", "yield_observed"] <- 3e11
plotYieldObservedVsModel(params)
```

Nothing breaks if your model keeps `yield_observed` among the species
parameters: a species that has no observation in the gear parameters takes its
value from there. What changes is that a model following mizer's own advice —
`given_species_params<-()` has been telling you to use `gear_params()<-` —
now works, where before the plot stopped with "You have not provided values for
the column 'yield_observed'".

### `matchYields()` and `calibrateYield()` have been removed

Both were deprecated in mizer 2.6.0 and nobody reported a use for them. They
adjusted the *abundance* of a species to move its yield, which is the wrong
lever: the yield is what the model predicts from the abundance and the fishing.
Replace `matchYields()` with `mizerExperimental::matchYield()`, which adjusts
the catchability instead:

```r
# Old
params <- calibrateYield(params)
params <- matchYields(params)
# New
params <- mizerExperimental::matchYield(params)
```

`calibrateYield()` has no replacement. It rescaled the whole model so that the
total yield summed over all species matched the total observation. If you were
using it to set the scale of your model, use `calibrateBiomass()` with observed
biomasses, or `scaleModel()` with a factor of your own choosing.

### The cheatsheet articles are now called guides

The topic articles that used to be called cheatsheets are called guides. A
cheatsheet reminds you of something you already know; these articles assume no
prior knowledge, so the name was wrong. Each article is now named after the
agent skill it is generated from, so that a topic has one name rather than
three, and its title is that skill's own heading:

| Old article | New article | New title |
|---|---|---|
| `cheatsheet-size-spectrum-dynamics` | `guide-understand-size-spectrum-dynamics` | Guide: Understanding size-spectrum dynamics |
| `cheatsheet-model-setup` | `guide-build-model` | Guide: Building a mizer model |
| `cheatsheet-calibration` | `guide-calibrate-model` | Guide: Reaching steady state and calibrating |
| `cheatsheet-changing-parameters` | `guide-change-parameters` | Guide: Changing model parameters |
| `cheatsheet-fishing` | `guide-set-up-fishing` | Guide: Setting up fishing |
| `cheatsheet-running-simulations` | `guide-run-simulation` | Guide: Running a mizer simulation |
| `cheatsheet-analysis-and-plotting` | `guide-analyse-and-plot` | Guide: Analysing and plotting mizer results |
| `cheatsheet-stability` | `guide-analyse-stability` | Guide: Analysing dynamic stability |
| `cheatsheet-extending-mizer` | `guide-extend-mizer` | Guide: Extending mizer |

"Using mizer extension packages" and "Creating a mizer extension package" are
now generated from skills too, so they are named after those skills like the
rest:

| Old article | New article | New title |
|---|---|---|
| `using-extension-packages` | `guide-use-extension-packages` | Guide: Using mizer extension packages |
| `creating-extension-packages` | `guide-create-extension-package` | Guide: Creating a mizer extension package |

The packaging article became a skill so that an agent helping you package an
extension can find it; it was previously the only extension document that was
not generated from one. Everything in the `extend-mizer` skill that only matters
once you share an extension moved into it at the same time, so the articles split
along that line: the mechanisms for changing mizer's dynamics in
**guide-extend-mizer**, and everything about turning that into a package other
people can install in **guide-create-extension-package**.

Its advice on marker classes was also corrected. It still told you to define
them with `setClass("myExtension", contains = "MizerParams")`, which mizer
3.2 made unnecessary and which actively prevents your package from being chained
with another, because a sealed class cannot be re-parented into the chain. Let
mizer create the classes; see the `create-extension-package` skill.

"Extending mizer" and "Guide: Extending mizer" were two articles on one topic,
the guide a short companion to the article. They are now a single guide,
generated from the `extend-mizer` skill, holding both the article's worked
examples and the guide's rules on quadrature schemes and discontinuous rates:

| Old article | New article | New title |
|---|---|---|
| `extending-mizer` | `guide-extend-mizer` | Guide: Extending mizer |

On the website the old addresses redirect, so a bookmark or a link in your own
writing still works. In R the old name does not resolve, because a vignette is
looked up by exactly its file name:

```r
# Old
vignette("cheatsheet-fishing")
# New
vignette("guide-set-up-fishing")
```

The `build-multispecies-model` skill was renamed to **build-model** in the same
pass: it covers `newTraitParams()`, `newCommunityParams()` and
`newSingleSpeciesParams()` as well, so its name claimed a narrower scope than it
has. If you install mizer's skills with `mizerAgents::setup_mizer_agent()`,
re-run it to pick up the new name.

### `projectToSteady()` ignores initial transients

To decide whether a simulation has settled onto a limit cycle,
`projectToSteady()` calculates the autocorrelation of a fine-resolution biomass
series. Previously it used the entire history from the start of the run. A large
initial transient could therefore dominate the autocorrelation and obscure a
cycle that had settled more recently.
In mizer 3.3, the autocorrelation step uses only the second half of the series
(or the most recent 20 samples if the series is shorter). This allows it to
ignore the initial transient. A cycle will now be found earlier (because the
check does not wait for the long-settled cycle to outweigh the transient), and
some cycles that were previously missed entirely will now be correctly reported.

