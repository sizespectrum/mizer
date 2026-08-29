## Upgrading from mizer 3.3 to 3.4

This release is a bug-fix release. Nothing moves unless your code calls
`scaleModel()`, or one of the calibration functions built on it, and then does
something that makes mizer recalculate the species parameters.

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
