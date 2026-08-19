## Upgrading from mizer 3.0 to 3.1

Version 3.1 leaves default results unchanged from 3.0 unless you opt in to the
new experimental second-order-in-size scheme. The changes below can still affect
existing code in specific situations.

### Maximum-size species parameters clarified

The maximum-size parameters have been given clearer, separate roles (#325):

- `w_inf`, the von Bertalanffy asymptotic size, is now the primary maximum-size
  parameter and is used as the default for `w_repro_max` (the size at which a
  mature individual invests all its energy in reproduction) and for `w_mat`.
- `w_max` is now purely a computational boundary — it sets the size grid and the
  plot range — and defaults to `1.5 * w_inf`.
- The default external mortality parameter `z0` is now computed from `w_inf`
  rather than `w_max`, so the computational boundary `w_max` no longer feeds into
  any model parameter.

Existing models and scripts are unaffected: if `w_inf` is not supplied it is
taken from `w_repro_max` or `w_max`, so old objects behave as before. However,
**new models built from the defaults may differ from 3.0.0**. If you build models
from scratch, check that `w_inf`, `w_max` and `w_repro_max` mean what you intend.

### `getTrophicLevel()` gives the resource a size-dependent trophic level

`getTrophicLevel()` and `getTrophicLevelBySpecies()` now assign the resource a
size-dependent trophic level,
$T_R(w) = \max(1,\, 1 + \log(w / w_R) / \log(\beta_R))$, instead of treating the
resource as trophic level 0. The new `w_R` and `beta_R` arguments control this.
Trophic levels computed with these functions will therefore be higher than
before. Set the arguments explicitly if you need to reproduce old numbers.

### Bug fixes that change results

Several fixes correct earlier behaviour and so change output:

- **`summary()` of a `MizerSim`** now reports the fishing effort that was used
  during the simulation, rather than the model's `initial_effort`. Gears whose
  effort varied over time show the mean, flagged with a note giving the range.
  The printed summary therefore differs for simulations run with time-varying
  effort.
- **`MizerSim` method for `plotDiet()`** introduced in version 3.0 simply
  plotted the diet at the initial time of the simulation. Now `plotDiet()` for
  a `MizerSim` accepts a `time_range` argument. The diet is now computed from
  the *simulated* abundances at the requested times, defaulting to the final
  saved step, rather than the initial one (#357).
- **Other components and `t_save`.** `project()` was advancing the abundances of
  other components (set via `setComponent()`) only once per *saved* time step
  instead of once per `dt` step. They are now integrated with the same `dt` as
  the consumer and resource spectra, so results for models with other components
  no longer depend on `t_save`.
- **Time-varying effort in `getRDI()`, `getRDD()`, `getFlux()`.** On a
  `MizerSim` object these now use the simulated time-varying effort rather than
  the initial effort, so they change for simulations with varying effort (#370).
- **`plotCDF()` / `plotlyCDF()` bin placement.** Each cumulative value is now
  plotted at its bin's *upper* edge, correcting a one-bin offset.
  The curves shift by one bin compared with 3.0 (#383).
- **`distanceMaxRelRDI()`.** Now returns `Inf` instead of `NaN` when a previous
  RDI is zero, so `projectToSteady()` no longer mistakes a `NaN` distance for
  convergence. Convergence behaviour can therefore differ in edge cases.

### Second-order methods advance the resource at the midpoint

If you use `project()` with `method = "predictor_corrector"` (or the new
`method = "tr_bdf2"`), the resource and the other components are now advanced
with midpoint rates rather than the start-of-step value, so that they reach the
same second-order accuracy in time as the consumer spectra. Results from these
methods therefore differ slightly from 3.0. The default `method = "euler"` and
the steady states are unchanged.

### Opting in to the second-order-in-size scheme

3.1 adds an optional, experimental second-order-accurate finite-volume scheme in
the size variable, controlled by the new `second_order_w` slot. It is **off by
default**, so default results are byte-identical to 3.0. If you switch it on
(via `second_order_w()<-` or the `second_order_w` argument of the `new...Params()`
constructors), size-integrated diagnostics and the resource spectrum shift by
$O(\Delta w)$, so a calibrated model may need recalibrating. See `?second_order_w`
and the "Numerical Details" vignette.

