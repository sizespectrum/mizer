## Upgrading from mizer 3.4 to 3.5

A model now remembers the edition of default choices it was created with. This
changes nothing for code that never calls `defaults_edition()`. Code that does
switch editions finds that the switch now applies only to models created
afterwards. Models saved with mizer 3.4 are upgraded once more when they are
first used. `plotYieldVsF()` accepts the `tol` argument it had in
mizerExperimental again, with a deprecation warning. On a model whose encounter
rate has been changed, `getDiet()` and `getTrophicLevel()` now follow the
change and say how they split it between the food sources. Projections with
high-order flux now report large growth Courant numbers; they retain the
previous fixed-step calculation unless subcycling is requested.

### A model keeps the defaults edition it was created with

`defaults_edition()` chooses the edition of the defaults that mizer fills in for
parameters you do not supply. Until now the choice was global: whenever mizer
filled in a default, for any model, it used the edition active at that moment.
An existing model therefore picked up the new edition's defaults as soon as
anything made mizer fill one in again, for example a change to a species
parameter, which recalculates every species parameter that mizer calculated.

Now the edition that is active when a model is created is recorded in the
model, and mizer always uses the model's own edition for it. Results move only
if your code activates one edition and then works with a model of another, most
often a model created or loaded before the `defaults_edition()` call, or a
model that comes with a package, such as `NS_params`, which is of edition 1.
Such a model now gets the defaults of its own edition when

- mizer recalculates a `gamma` that it calculated, where edition 2 takes
  `interaction_resource` into account and edition 1 does not,
- `setReproduction()` recalculates `psi`, which edition 1 sets to 1 above
  `w_repro_max`,
- you assign gear parameters without a `catchability`,
- you give `project()` or `initial_effort<-()` an effort that leaves out some
  gears or contains `NA`,
- `addSpecies()` adds species and gears to it,
- `get_initial_n()` calculates initial abundances for it.

See the edition of a model with `defaults_edition(params)`. To give an existing
model the defaults of a different edition, create it again with that edition
active.

A model saved with an earlier version of mizer does not record its edition, and
mizer cannot tell from the model which edition was active when it was created.
It treats every such model as edition 1, which was the default. If you created
yours with edition 2 active, say so once after loading it:

```r
defaults_edition(params) <- 2
```

This recalculates nothing; it only tells mizer which defaults to use for the
model from now on.

Recording the edition is a change to the structure of a model, so a model saved
with mizer 3.4 is upgraded when it is first used, with the usual warning that it
was created with an earlier version of mizer. `readParams()` and `readSim()` do
the upgrade quietly. Reassign the result, `params <- validParams(params)`, to
stop the warning, and save the model again to keep the upgraded version.

### `plotYieldVsF()` accepts `tol` again

When `plotYieldVsF()` moved to mizer from mizerExperimental in mizer 3.3, its
`tol` argument became `distance_tol`, and in mizer 3.3 and 3.4 a call that
still passed `tol` stopped with an error. It now runs, with a deprecation
warning, and uses the value as `distance_tol`, which is what it always meant.
Rename the argument to stop the warning:

```r
# Old
plotYieldVsF(params, "Cod", F_range = seq(0.2, 0.8, 0.02), tol = 1e-4)
# New
plotYieldVsF(params, "Cod", F_range = seq(0.2, 0.8, 0.02), distance_tol = 1e-4)
```

Giving both `tol` and `distance_tol` is an error.

### `getDiet()` and `getTrophicLevel()` follow a changed encounter rate

On a model whose encounter rate has been changed, by a function you registered
with `setRateFunction("Encounter", ...)` or by an extension package,
`getDiet()`, `plotDiet()` and `getTrophicLevel()` now give a message that starts

```
The encounter rate has been changed by a custom encounter rate function or by an
extension, and `getEncounterByPrey()` does not say how the change is shared
between the food sources.
```

These functions used to ignore the change, so the diet did not add up to what
the predators consumed and the trophic levels could be far off. They now take
the total from the model's own encounter rate. What they cannot know is how
your change is shared between the prey species, the resource and the other food
sources, so they scale all of them in proportion. That is exact if your
function multiplies the whole encounter rate by a factor, as a seasonal or
temperature multiplier does, and the message is then nothing to act on. If your
function adds food of its own, the extra food is spread over the existing food
sources in the diet, or shown as "External" at sizes at which they provide
nothing. Add it as a component with an `encounter_fun` instead, see
`setComponent()`, and `getDiet()` shows it as a food source of its own. If the
message comes from an extension package, a newer version of that package can
supply the exact split.

To silence the message, set `options(mizer_info_level = 0)` or wrap the call in
`suppressMessages()`.

### Projection reports large growth Courant numbers

With van Leer or centred flux, `project()` and projection-based settling can
now warn that the growth Courant number exceeds the advisory target. Existing
scripts may see this warning even if they previously ran without a report.
The default timestep calculation has not changed. The warning identifies a
risk of alternating cell-scale numerical residuals, which small changes in
total biomass can conceal.

Inspect the current state with `getCourant(params, dt = dt)`, reduce `dt`, or
enable coupled subcycling:

```r
sim <- project(params, t_max = 20, dt = 1 / 52, method = "tr_bdf2",
               subcycle = TRUE, max_courant = 1)
```

Subcycling advances fish, resources and other components together with
smaller internal steps. It can substantially increase runtime. The stored
settings carry over when continuing the simulation with `project(sim)`.
The target is growth-based guidance, not a guarantee of stability or
positivity for all dynamics. Check results under timestep refinement.

If the timestep has been checked and the report is unwanted, use
`info_level = 0` for the call, or set `options(mizer_info_level = 0)`.


### Plot limits only change the visible window

`wlim`, `llim`, `tlim` and `ylim` now zoom the plot without filtering its data.
This applies to array plots, named plots, comparisons and relative plots.
`return_data = TRUE` returns the same data with or without these limits, and
community totals and cumulative distributions no longer depend on the window.
Lines and areas can cross a window boundary without losing their neighbouring
vertices. Animations still use `tlim` to choose their displayed frames.

Code that used plot limits to extract a subset of the returned data should
subset that data explicitly. For biomass or abundance integrated over a size
range, continue to use the getter's `min_w`/`max_w` or `min_l`/`max_l` arguments.

For `plotCDF()`, `plotCDF2()` and their `plotly…()` counterparts, select the
integration range with the new scalar `min_w`, `max_w`, `min_l` and `max_l`
arguments. Length bounds override the corresponding weight bounds and are
converted with each species' own length-weight parameters, including the
resource's equivalent-sphere convention. Partial boundary bins contribute only
the part inside the interval, under both quadrature schemes.

```r
# Previously used the visible size range to select a distribution:
plotCDF(params, wlim = c(10, 1000))
# Select the distribution explicitly, and optionally zoom to the same range:
plotCDF(params, min_w = 10, max_w = 1000, wlim = c(10, 1000))
```

Without integration bounds the CDF covers the complete model grid, including
resource sizes previously hidden by the spectrum plot's default window.
Totals include all species and the resource regardless of which are displayed.
A completed cumulative contribution stays constant beyond its last bin, so a
community CDF on a length axis does not fall when one species' grid ends.
