# Project to steady state

**\[experimental\]**

Run the full dynamics, as in
[`project()`](https://sizespectrum.org/mizer/reference/project.md), but
stop once the change has slowed down sufficiently, in the sense that the
distance between states at successive time steps is less than `tol`. You
determine how the distance is calculated.

Because reproduction stays dynamic here, the run can only ever end up on
an attractor of the dynamics, and that need not be a fixed point: it may
stop on a limit cycle, on a species going extinct, or simply at `t_max`.
So the state it leaves behind is not necessarily a steady state; see the
section below on how to check.

## Usage

``` r
projectToSteady(
  params,
  effort = params@initial_effort,
  distance_func = distanceSSLogN,
  t_per = 1.5,
  t_max = 100,
  dt = 0.1,
  t_save = dt,
  tol = 0.1 * t_per,
  amplitude_tol = 0.01,
  amp_rel_tol = 0.1,
  extinction_threshold = 1e-06,
  return_sim = FALSE,
  progress_bar = TRUE,
  info_level = default_info_level(),
  method = c("euler", "predictor_corrector", "tr_bdf2"),
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- effort:

  The fishing effort to be used throughout the simulation. This is
  validated by
  [`validEffortVector()`](https://sizespectrum.org/mizer/reference/validEffortVector.md)
  and can therefore be `NULL`, a single numeric value used for all
  gears, an unnamed numeric vector with one entry per gear, or a named
  numeric vector for some or all gears.

- distance_func:

  A function that will be called after every `t_per` years with both the
  previous and the new state and that should return a number that in
  some sense measures the distance between the states. By default this
  uses the function
  [`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md)
  that you can use as a model for your own distance function.

- t_per:

  The simulation is broken up into shorter runs of `t_per` years, after
  each of which we check for convergence. Default value is 1.5. This
  should be chosen as an odd multiple of the timestep `dt` in order to
  be able to detect period 2 cycles.

- t_max:

  The maximum number of years to run the simulation. Default is 100.

- dt:

  The time step to use in
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- t_save:

  The interval at which a cheap per-species biomass summary is recorded
  for limit-cycle detection. Must be a positive multiple of `dt` and a
  divisor of `t_per`. Smaller values resolve the cycle period more
  finely at a small extra cost. Default is `dt`.

- tol:

  The run stops when the number returned by `distance_func` for two
  states `t_per` years apart drops below `tol`. Its meaning therefore
  depends on the distance function you supply.

- amplitude_tol:

  **\[experimental\]** The minimum relative biomass amplitude for a
  persistent oscillation to be reported as a limit cycle rather than
  treated as an (effectively steady) fixed point. This is a fraction of
  mean biomass and is kept separate from `tol` (which measures
  convergence to a fixed point on a different scale). Default `0.01`.

- amp_rel_tol:

  **\[experimental\]** Maximum relative change of amplitude between
  successive periods for the cycle to count as settled. Default `0.01`.

- extinction_threshold:

  **\[experimental\]** A species is treated as going extinct, stopping
  the run, once its reproduction rate (RDD) falls below this fraction of
  its value at the start of the run. For example the default `1e-6`
  treats a species as extinct once its reproduction has collapsed to a
  millionth of its initial level. Because it is relative to the initial
  reproduction, a species that starts with zero reproduction is flagged
  immediately, and (in
  [`steady()`](https://sizespectrum.org/mizer/reference/steady.md),
  where reproduction is held constant) a healthy species is never
  flagged.

- return_sim:

  If TRUE, the function returns the MizerSim object holding the result
  of the simulation run, saved at intervals of `t_per`. If FALSE
  (default) the function returns a MizerParams object with the "initial"
  slots set to the steady state.

- progress_bar:

  A shiny progress object to implement a progress bar in a shiny app.
  Default FALSE.

- info_level:

  Controls the amount of information messages that are shown. Higher
  levels lead to more messages, `info_level = 0` gives silence. The
  default is taken from the `mizer_info_level` option, see
  [`default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.md).

- method:

  The numerical method to use for the consumer density update. See
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- ...:

  Further arguments will be passed on to your distance function.

## Value

If `return_sim = FALSE`, a `MizerParams` object with the initial state
replaced by the final state found by the steady-state search. If
`return_sim = TRUE`, a `MizerSim` object containing the intermediate
states saved every `t_per` years.

In either case the returned object carries an attribute `"convergence"`
describing the solution the run settled on, a named list with entries:

- `type`:

  One of `"steady"` (a stable fixed point), `"cycle"` (a limit cycle),
  `"not_converged"` (still changing at `t_max`) or `"extinction"` (a
  species died out).

- `converged`:

  Logical, `TRUE` for `"steady"` and `"cycle"`.

- `distance`:

  The final value returned by `distance_func`.

- `residual`:

  The largest per-capita rate of change, in 1/year, at the state that
  was reached, as returned by
  [`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md).
  Unlike `distance`, which compares two states `t_per` apart on whatever
  scale the distance function uses, this measures how far the state
  actually is from being a fixed point.

- `years`:

  The number of years simulated.

- `period`:

  For a limit cycle, its period in years; otherwise `NA`.

- `amplitude`:

  For a limit cycle, the largest per-species relative peak-to-trough
  biomass amplitude; otherwise `NA`.

This mirrors how
[`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md)
attaches an `"stability"` attribute.

## How the run is organised

The simulation is not run in one go but is broken into blocks of `t_per`
years. Within a block the dynamics are advanced with time step `dt`
exactly as in
[`project()`](https://sizespectrum.org/mizer/reference/project.md). At
the end of each block the function decides whether to stop, so `t_per`
sets how often the stopping criteria are evaluated and also the interval
over which change is measured. The run ends at the latest after `t_max`
years, i.e. after `floor(t_max / t_per)` blocks.

Independently of the blocks, a cheap scalar summary of the state (the
biomass of each species) is recorded every `t_save` years. This finely
resolved series is what the limit-cycle detection works on, so that a
cycle can be found and its period measured even when that period bears
no simple relation to `t_per`.

At the end of each block the following checks are made, in this order.

### 1. Extinction

If the reproduction rate (RDD) of any species has fallen below
`extinction_threshold` times its value at the start of the run, or has
become `NA`, that species is deemed to be on its way to extinction. A
warning naming the affected species is issued and the run stops with
`type = "extinction"`. Because the criterion is relative to the initial
reproduction, a species that starts with zero reproduction is flagged
immediately, whereas in
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md), where
reproduction is held constant, a healthy species is never flagged.

### 2. Convergence to a fixed point

`distance_func` is called with the state at the end of the previous
block and the state at the end of the current block, i.e. with two
states `t_per` years apart. If the number it returns is less than `tol`,
the run stops with `type = "steady"`. What "distance" means is entirely
up to that function: the default
[`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md)
uses the sum of squared changes in log abundance, while
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) instead
passes
[`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md),
which uses the largest relative change in egg production.

Note that this test only compares two states one `t_per` apart; it
cannot by itself tell a fixed point from a cycle whose period divides
`t_per`. This is why `t_per` should be chosen as an *odd* multiple of
`dt`: a period-2 cycle (period `2 * dt`), which is the most common
numerical oscillation, would otherwise be sampled at the same phase in
every block and would look perfectly converged.

### 3. Limit cycle

If the state is still changing more than `tol`, the recorded biomass
series is examined to see whether the run has instead settled onto a
limit cycle. If it has, the run stops with `type = "cycle"` and the
period and amplitude of the cycle are reported.

If none of the three checks fires before `t_max` is reached, the run
stops with `type = "not_converged"`. In every case the outcome is
recorded in the `"convergence"` attribute of the returned object,
described under *Value* below.

## How a limit cycle is detected

The detection uses the community-total biomass, on a log scale and with
its mean removed, as a scalar signal, sampled every `t_save` years. At
least 20 samples are needed before any cycle can be reported.

1.  **Candidate period.** The autocorrelation function of the signal is
    computed up to a lag of half the length of the series, and the first
    local maximum with an autocorrelation above `0.5` is taken as the
    candidate period. If there is no such peak, or the peak is at a lag
    of one sample, no cycle is reported.

2.  **Enough history.** The series must cover at least three full
    candidate periods. Otherwise the check is deferred to a later block,
    when more history has accumulated.

3.  **Amplitude.** For each of the last three period-long windows, the
    amplitude is measured as the largest over species of the relative
    peak-to-trough biomass range `(max - min) / mean`. The amplitude in
    the most recent window must exceed `amplitude_tol`; a smaller
    oscillation is considered negligible and the state is left to be
    treated as a fixed point.

4.  **Settled.** The amplitudes of the three successive windows must
    agree with each other to within `amp_rel_tol`, and the most recent
    amplitude must not be smaller than the oldest by more than
    `amp_rel_tol`.

The last condition is what distinguishes a genuine limit cycle from a
slowly decaying spiral towards a stable fixed point: the spiral loses
amplitude from one period to the next, the cycle does not. The
distinction is necessarily imperfect when the decay is extremely slow,
because over any finite run such a spiral is indistinguishable from a
cycle. If you need a definitive answer, use
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
on the fixed point found by
[`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md),
which works out the spectral radius of the linearised dynamics instead
of watching a trajectory.

The reported `period` is a multiple of `t_save`, so it is only resolved
to that accuracy; reduce `t_save` if you need the period more precisely.

## What you get back may not be a steady state

The stopping criterion is a proxy. It says that two states `t_per` years
apart differ by less than `tol` on whatever scale the criterion is
measured on; it does not say that the state reached is a fixed point.
There are four ways the returned object can fail to be one:

- the run converged on its own scale while the biomasses are still
  visibly drifting, because `tol` was too loose;

- the run reached `t_max` without converging (`type = "not_converged"`);

- the run settled on a limit cycle (`type = "cycle"`), in which case the
  state stored is one point on that cycle;

- the run stopped because a species was going extinct
  (`type = "extinction"`).

So treat the result as a claim to be checked rather than as a guarantee:

    attr(params, "convergence")$type      # "steady", "cycle", "not_converged", "extinction"
    attr(params, "convergence")$residual  # largest biomass drift, in 1/year
    isSteady(params)                      # TRUE if within tolerance
    summary(params)                       # includes the biomass-drift verdict
    plot(getSteadyResidual(params))       # which species, and at which sizes

The first four say *whether* the model is steady, the last says *where*
it is not, which is the one to reach for when it is not: a model that is
off steady state is usually off in one species or one part of the size
range, and the plot names it. See
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md)
for why the verdict is phrased in terms of biomass drift rather than the
largest per-capita rate.

The messages this function prints say the same thing — a converged run
whose biomasses are still moving reports the drift and adds "Reduce
`tol` to converge further." — but they are suppressed by
`info_level = 0`, so in a script the `"convergence"` attribute is the
reliable check.

Finally, a genuine fixed point need not be a *stable* one. Use
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
to find out, and
[`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md)
to converge onto a fixed point that the dynamics themselves would run
away from.

## See also

[`steady()`](https://sizespectrum.org/mizer/reference/steady.md),
[`isSteady()`](https://sizespectrum.org/mizer/reference/isSteady.md),
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md),
[`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md),
[`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md),
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
