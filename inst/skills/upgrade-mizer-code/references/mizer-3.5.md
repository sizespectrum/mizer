## Upgrading from mizer 3.4 to 3.5

A model now remembers the edition of default choices it was created with. This
changes nothing for code that never calls `defaults_edition()`. Code that does
switch editions finds that the switch now applies only to models created
afterwards. Models saved with mizer 3.4 are upgraded once more when they are
first used. `plotYieldVsF()` accepts the `tol` argument it had in
mizerExperimental again, with a deprecation warning.

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
