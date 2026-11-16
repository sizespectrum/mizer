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

### Component updates retain omitted settings

`setComponent()` can update an existing component without repeating its initial
state or dynamics function. Omitted component parameters, callbacks, colour and
line type are retained. Earlier versions reset omitted plotting styles to grey
and solid on each update; supply `colour = "grey", linetype = "solid"` explicitly
if that reset is wanted.

The new `definition` argument supplies a reusable named list of callbacks and
plotting styles. Passing a definition replaces the old definition before explicit
arguments override its entries. Omitted optional callbacks are then removed,
and omitted styles take the defaults; state and component parameters remain
unchanged unless supplied separately. Explicit `NULL` removes an optional
encounter or mortality callback.

`getComponent()` now returns `definition`, `rate_fun`, `convergence_fun` and
`role` entries alongside its original five entries. Code that selects entries
by name continues to work;
code that relies on an exact five-element result should select the entries it
needs. The assembled definition is saved in the model, so changing the original
list does not silently update models that already used it.

### Component correctors use complete midpoint rates

Under `method = "predictor_corrector"` or `"tr_bdf2"` (also called
`"second_order"`), custom resource and component steps now receive midpoint
values for every standard numeric rate. Previously only growth, mortality,
diffusion, recruitment and resource mortality were averaged; other entries such
as `pred_rate` and `f_mort` retained their start values. Results may therefore
change for models whose custom dynamics use those entries. Euler stepping is
unchanged.

Components can also supply equations for mizer to integrate with Euler or Heun.
The fish solver does not determine a custom component's numerical order or
stability. Inspect `sim@sim_params$component_integration` for the recorded method
and declared order of each component; an unknown custom order is shown as `NA`.

### Steady-state checks cover the requested model scope

`isSteady()` now has three possible answers. `TRUE` means every required measure
passed; `FALSE` means at least one failed; `NA` means no failure was established
but required support was missing. Use `isTRUE(isSteady(params))` when code must
proceed only after a positive verdict. Inspect `getSteadyAssessment(params)` for
the measures, tolerances and missing capabilities. Callback errors now propagate
instead of being silently turned into a negative or unavailable verdict.

The default built-in criterion integrates the magnitudes of size-class changes,
using the model's quadrature and documented absolute and relative tolerances.
Changes at different sizes can no longer cancel. A model with little net biomass
change may therefore fail this stronger check. Request `criterion = "biomass"`
explicitly to use net biomass drift. `getSteadyResidual()` keeps its detailed
residual-array meaning.

Steady-state checks now include components by default. A component without an
assessment leaves coverage incomplete, even if projection works. Read-only
assessment reports this; a solver requiring that support stops before starting.
`project(check_steady = TRUE)` warns about incomplete coverage, and `summary()`
reports the assessment and dimensionless tolerance ratios instead of using one
biomass number as a whole-model verdict.

To request a subsystem explicitly, use `scope = "conditional"` and `fixed`, for
example `fixed = "component:pool"`. Conditional solver results retain their
scope, and a final restored model is reassessed. The current Newton solver cannot
evolve arbitrary dynamic components; it requires those states to be explicitly
fixed. The built-in `residual` convergence field remains a biomass diagnostic;
read `assessment` for the verdict behind `attractor`.

### Tuning uses declared balancing strategies

`tuneSteadyState()` now uses the default strategy stored for each process.
For the usual Beverton-Holt model this still holds resource abundance and
recruitment fixed and preserves reproduction level when balancing afterwards.
Other laws use their own defaults. Supply `strategies` to select alternatives,
or `joint` for a registered strategy covering coupled processes. An explicit
`preserve` request must be supported by the selected reproduction law.

Missing stationary tuning support stops before the search. A model that can
project is not necessarily one that can be tuned; update the extension package
when it reports missing support. Mizer checks that balancing preserves the
target state and every effective rate entering its equations. It no longer
accepts a different set of rates simply because their net residual cancels,
and does not retry with a generic optimiser or repair the result by solving
again. `search_assessment` describes the constrained search; `assessment`
describes the restored model.

`balanceRDD()` retains the current reproduction law. With no target it preserves
current effective recruitment; explicit `rdd_target` values are useful for
calibration. Matching growth, biomass or abundance no longer implicitly changes
the law to Beverton-Holt, and partial selections leave unselected reproduction
parameters alone. Built-in overcompensatory laws use the ascending branch by
default; request the descending branch explicitly when that is intended.

`reproduction_level()` is now law-specific. It is recruitment divided by the
law's finite maximum when one exists. No density dependence and unbounded
Shepherd curves have level zero; constant positive recruitment has level one,
and constant zero recruitment has an undefined (`NA`) level.

### Beverton-Holt switching respects explicit controls

`setBevertonHolt()` is now an explicit switch to that reproduction law, retaining
the model's effective recruitment at its current state. On a model away from
steady state this can differ from the recruitment required to hold egg abundance
constant. To request that separate target, use
`balanceRDD(params, rdd_target = getRequiredRDD(params))` after selecting the law.

An explicitly supplied reproductive efficiency, capacity or reproduction level
must be feasible and is no longer silently clamped or changed. In particular,
positive recruitment cannot be produced by an absent spawning population, and a
positive reproduction level cannot describe a species with zero recruitment.
Choose a feasible control or change the underlying model state deliberately.

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
