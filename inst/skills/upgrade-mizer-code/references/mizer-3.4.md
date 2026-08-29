## Upgrading from mizer 3.3 to 3.4

Nothing moves unless your code calls `scaleModel()`, or one of the calibration
functions built on it, and then does something that makes mizer recalculate the
species parameters, or unless you load a saved model whose species parameters
had been edited by writing into the `species_params` slot directly.

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
when you have just scaled those arrays by hand yourself.
