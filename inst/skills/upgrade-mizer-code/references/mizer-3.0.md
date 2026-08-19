## Upgrading from mizer 2.5.4 to 3.0

Version 3.0 is a large release. Most new capabilities are additive and off by
default, but there are several renamed arguments, deprecations and behavioural
changes that can affect existing code.

### Renamed arguments and changed defaults (breaking changes)

- **First argument of `plotBiomass()`, `plotYield()`, `plotYieldGear()`** (and
  their `MizerSim` methods and `plotly*` wrappers) is renamed from `sim` to
  `object`, for consistency with the other plot generics. Calls that passed the
  simulation by name, `plotBiomass(sim = my_sim)`, must become
  `plotBiomass(object = my_sim)`. Positional calls are unaffected.
- **`plotBiomassObservedVsModel()` / `plotlyBiomassObservedVsModel()`** now
  default to `ratio = FALSE` for all object types. Calls that relied on the
  previous ratio plot must set `ratio = TRUE` explicitly.
- **`plotDiet()` no longer accepts a `time_range` argument.** Remove it from your
  calls. (In 3.1 a `time_range` argument returns for the `MizerSim` method — see
  above.)
- **Dimnames of `getMort()` and `getPredRate()`** arrays are now `sp` and `w`
  (matching `getFMort()` and the other rate getters). Code that referred to the
  old dimnames by name must be updated.

### Rate getters return classed array objects

Functions that return arrays of the form (species × size), (time × species) or
(time × species × size) now attach extra attributes and an S3 class
(`ArraySpeciesBySize`, `ArrayTimeBySpecies` or `ArrayTimeBySpeciesBySize`). The
numeric values and ordinary matrix behaviour (arithmetic, subsetting) are
unchanged, but the extra class and attributes mean that a strict comparison such
as `identical(getMort(params), old_value)` can now report a difference where the
numbers agree. Use `unclass()`, or compare with `all.equal()` on the values, if
you need to ignore the class. These objects also carry `print()`, `summary()`,
`plot()` and `as.data.frame()` methods, so printing them looks different from a
bare matrix.

### `setInitialValues()` is deprecated

`setInitialValues()` is deprecated. Replace

```r
params <- setInitialValues(params, sim)
```

with

```r
params <- finalParams(sim)
```

or, when averaging over a time range, with
`getParams(sim, time_range, geometric_mean)`. This reflects a shift in
interpretation: a `MizerParams` object now represents not just the model
specification but also its current state (the abundances), which can be
extracted from a simulation with `getParams()`, `finalParams()` and
`initialParams()`.

### Growth can no longer be negative

Growth is now forced to be non-negative, preventing unphysical shrinkage. In any
model where the energy available for growth used to go negative (for example a
strongly food-limited large individual), growth is now clamped at zero instead,
so projected size spectra can differ from 2.5.4. No warning is issued when growth
stops at or after the maturity size.

### `project()` timing and effort handling

- **Inherited `dt` and `method`.** When `project()` is called on an existing
  `MizerSim` object, `dt` and `method` now default to the values stored in the
  simulation's new `sim_params` slot. If you pass values that differ from the
  stored ones, a warning is issued. To use different settings deliberately, pass
  them explicitly and expect the warning.
- **`t_max` / `t_save` with an effort array.** These arguments are now respected
  even when an effort array is supplied (#231). With `t_max` the simulation
  extends beyond the times in the effort array using the last known effort; with
  `t_save` the save frequency is controlled independently, interpolating effort
  as needed. Simulations that previously derived their length or save times
  solely from the effort array may now produce a different set of saved steps.
- **State at `t_max` always saved.** `project()` now warns when `t_max` is not a
  multiple of `t_save` and ensures the state at `t_max` is saved even if the
  final interval is shorter than `t_save` (#341). The returned simulation may
  therefore contain one extra saved time step compared with 3.0.

### `plot()` and `summary()` are now S3 methods

The `plot()` and `summary()` methods for `MizerParams`, `MizerSim` and the mizer
array classes are now registered as S3 methods rather than S4 methods, so
`plot()` and `summary()` stay plain S3 generics when mizer is loaded. This avoids
interfering with S4 dispatch in other packages, but code that relied on
`plot`/`summary` being S4 generics (for example via `selectMethod()` or
`getMethod()`) needs adjusting.

### Bug fixes that change results

- **`getMeanMaxWeight()`** now applies the species selector to the denominator as
  well, so its values change when a subset of species is selected.
- **`plotSpectra()` axis limits.** It no longer forces the y-axis lower limit to
  `1e-20` (it auto-scales to the data) and, when `resource = FALSE`, it uses
  `min(params@w)` rather than `min(params@w) / 100` as the default lower size
  limit. Plots therefore look different.
- **`getFMort()` on a `MizerSim`** was silently dropping the component names from
  `n_other`, breaking rate functions that access `n_other` by name (e.g.
  `n_other[["resource"]]`); it now preserves them.
- **`getFMort.MizerSim()`** now passes the time argument `t` to user-defined
  fishing-mortality functions, so a time-dependent fishing function now sees the
  correct time.

### Predation diffusion is available but off by default

3.0 adds a diffusion term to the growth dynamics, controlled by the new
`use_predation_diffusion` slot. It defaults to `FALSE`, preserving the behaviour
of earlier mizer, so existing models are unchanged unless you switch it on with
`use_predation_diffusion(params) <- TRUE`. Likewise the new species parameters
`z_ext`, `d`, `E_ext` and `D_ext` for external mortality, encounter and diffusion
all default to values that leave the model unchanged.
