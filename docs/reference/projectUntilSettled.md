# Project the dynamics until they settle

**\[experimental\]**

Run the full dynamics, as in
[`project()`](https://sizespectrum.org/mizer/reference/project.md), but
stop once the run has settled: either the change has slowed down
sufficiently, in the sense that the distance between states `t_check`
years apart is less than `distance_tol` and the state has stopped
drifting, or the run has been recognised as being on a limit cycle. You
determine how the distance is calculated.

Nothing is held fixed, so the run can only ever end up on an attractor
of the dynamics, and that need not be a fixed point: besides a limit
cycle it may stop on a species going extinct, or simply at `t_max`.
"Settled" is therefore the most this function claims; the state it
leaves behind is not necessarily a steady state. See the section below
on how to check, and use
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
if what you want is the steady state itself rather than the trajectory
leading to it.

## Usage

``` r
projectUntilSettled(
  params,
  effort = params@initial_effort,
  distance_func = distanceSSLogN,
  t_check = 15 * dt,
  t_max = 100,
  dt = 0.1,
  t_save = 1,
  distance_tol = 0.1 * t_check,
  residual_tol = steady_residual_tol(),
  amplitude_tol = 0.01,
  amp_rel_tol = 0.1,
  extinction_threshold = 1e-06,
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

  A function that will be called at every check with both the previous
  and the new state and that should return a number that in some sense
  measures the distance between the states. By default this uses the
  function
  [`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md)
  that you can use as a model for your own distance function.

- t_check:

  The interval in years at which the run pauses to check whether it has
  settled, and hence also the interval over which `distance_func`
  measures change. Must be a positive multiple of `dt`. The default
  `15 * dt` is an odd multiple of the time step, which is what lets a
  period-2 cycle be seen; you should rarely need to change it.

- t_max:

  The maximum number of years to run the simulation. Default is 100.

- dt:

  The time step to use in
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- t_save:

  The interval in years at which the state is stored in the returned
  `MizerSim`, as in
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).
  Must be a positive multiple of `dt`, but need bear no relation to
  `t_check`. Default is 1. The state the run settles on is always the
  final time point, even when the run stops between two saves.

- distance_tol:

  The run stops when the number returned by `distance_func` for two
  states `t_check` years apart drops below `distance_tol`, provided the
  drift criterion below is also met. Its meaning therefore depends on
  the distance function you supply. It was called `tol` before mizer
  3.3.

- residual_tol:

  **\[experimental\]** The largest relative rate of biomass change, in
  1/year, that a state may have and still be called a fixed point. This
  is the criterion of
  [`isSteady()`](https://sizespectrum.org/mizer/reference/isSteady.md)
  and is measured with
  [`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
  so it has the same meaning for every model and every distance
  function. Default `0.05`.

  It is a backstop against a distance function that has gone quiet while
  the model is still moving, not the main line of defence against a
  cycle: an oscillation of relative amplitude \\A\\ and period \\T\\
  drifts at up to \\2\pi A/T\\ per year, which for a small
  `amplitude_tol` can be below any sensible `residual_tol`. Cycles that
  small are caught by the cycle detection above, which is why that runs
  unconditionally.

- amplitude_tol:

  **\[experimental\]** The minimum relative biomass amplitude for a
  persistent oscillation to be reported as a limit cycle rather than
  treated as an (effectively steady) fixed point. This is a fraction of
  mean biomass and is kept separate from `distance_tol` (which measures
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
  [`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md),
  where reproduction is held constant) a healthy species is never
  flagged.

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

A `MizerSim` object containing the states saved every `t_save` years,
with the state the run settled on as its final time point. That last
interval is shorter than the others when the run stops between two
saves. Use
[`finalParams()`](https://sizespectrum.org/mizer/reference/getParams.md)
to extract that final state as a `MizerParams` object, or call
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
instead, which returns it directly.

The returned object carries an attribute `"convergence"` describing the
solution the run settled on, a named list with entries. The first three
answer three different questions and should not be read as one:

- `termination`:

  Why the run stopped: `"residual_tolerance"` (both convergence criteria
  were met), `"distance_tolerance"` (the distance function was satisfied
  but the state is still drifting — reachable with a loose
  `residual_tol`, and from the superseded
  [`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.md),
  which stops on the distance criterion alone), `"cycle_detected"` (a
  limit cycle), `"time_limit"` (still changing at `t_max`) or
  `"extinction"` (a species died out). The steady-state finders can also
  return `"solver_converged"` and `"solver_failed"` from the Newton
  solver.

- `converged`:

  Logical, `TRUE` when the run stopped on a criterion of its own rather
  than running out of time or losing a species. This is a statement
  about the numerics, not about the state.

- `attractor`:

  What the state that was reached actually is: `"fixed_point"` when the
  biomass drift is within `residual_tol`, `"limit_cycle"` when a cycle
  was detected, and `NA` when it is neither — a run stopped in
  mid-flight, or a species on its way out. This is the entry to test
  before treating a result as a steady state.

- `distance`:

  The final value returned by `distance_func`.

- `residual`:

  The largest absolute relative rate of biomass change, in 1/year, at
  the state that was reached. For each consumer species this is
  \\(dB_i/dt) / B_i\\, a biomass-weighted aggregate over its size
  classes; the resource is treated the same way and the other components
  contribute their own relative rates. Unlike `distance`, which compares
  two states `t_check` apart on whatever scale the distance function
  uses, this measures how far the state actually is from being a fixed
  point.

- `years`:

  The number of years simulated. `NA` for a direct solve.

- `period`:

  For a limit cycle, its period in years; otherwise `NA`.

- `amplitude`:

  For a limit cycle, the largest per-species relative peak-to-trough
  biomass amplitude; otherwise `NA`.

- `extinct`:

  Character vector naming any species that went extinct during the run,
  or `character(0)` if none.

## How the run is organised

The dynamics are advanced with time step `dt` exactly as in
[`project()`](https://sizespectrum.org/mizer/reference/project.md).
Every `t_check` years the function pauses to decide whether to stop, so
`t_check` sets how often the stopping criteria are evaluated and also
the interval over which change is measured. You should not normally need
to set it: it defaults to `15 * dt`, which is an odd multiple of the
time step for the reason given below. The run ends at the latest after
`t_max` years.

Independently of that, the state is stored in the returned `MizerSim`
every `t_save` years, exactly as in
[`project()`](https://sizespectrum.org/mizer/reference/project.md), and
a cheap scalar summary of the state (the biomass of each species) is
recorded after every time step. That finely resolved series is what the
limit-cycle detection works on, so that a cycle can be found and its
period measured even when that period bears no simple relation to
`t_check`. The three intervals are independent of each other; `t_check`
and `t_save` need only be multiples of `dt`.

At each check the following tests are made, in this order.

### 1. Extinction

If the reproduction rate (RDD) of any species has fallen below
`extinction_threshold` times its value at the start of the run, or has
become `NA`, that species is deemed to be on its way to extinction. A
warning naming the affected species is issued and the run stops with
`type = "extinction"`. Because the criterion is relative to the initial
reproduction, a species that starts with zero reproduction is flagged
immediately, whereas in
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md),
where reproduction is held constant, a healthy species is never flagged.

### 2. Limit cycle

The recorded biomass series is examined to see whether the run has
settled onto a limit cycle. If it has, the run stops with
`type = "cycle"` and the period and amplitude of the cycle are reported.

It is made at every check, whether or not the state looks converged by
the measure below. A cycle whose period divides `t_check` puts the two
states that the distance function compares at the same phase, so it
would otherwise be reported as a fixed point of zero width. The
detection works on the biomass series sampled at every time step
instead, which is blind to `t_check`.

### 3. Convergence to a fixed point

Two things have to hold for the run to stop on a fixed point.

First, `distance_func` is called with the state at the previous check
and the state at the current one, i.e. with two states `t_check` years
apart, and the number it returns must be less than `distance_tol`. What
"distance" means is entirely up to that function: the default
[`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md)
uses the sum of squared changes in log abundance, while
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
instead passes
[`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md),
which uses the largest relative change in egg production.

Second, the state actually reached must be a fixed point: the largest
relative rate of biomass change there, as measured by
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
must be at most `residual_tol`. The distance criterion on its own says
only that the state stopped moving on the scale of the distance
function, which is a different question — a distance function can be
insensitive to the very motion that is left. When the distance criterion
is met but the drift is not, the run carries on rather than declaring a
fixed point.

When both hold the run stops with `termination = "residual_tolerance"`.
That is deliberately not called `"steady"`: `residual_tol` is a working
tolerance rather than a proof, and the `residual` entry of the
`"convergence"` attribute reports the drift that was actually reached.

Even so, `t_check` should be an *odd* multiple of `dt`, which is why it
defaults to `15 * dt`. A period-2 cycle (period `2 * dt`), the most
common numerical oscillation, is otherwise sampled at the same phase at
every check, and its amplitude can sit below `amplitude_tol` where the
cycle detection deliberately ignores it.

If none of the three checks fires before `t_max` is reached, the run
stops with `type = "not_converged"`. In every case the outcome is
recorded in the `"convergence"` attribute of the returned object,
described under *Value* below.

## How a limit cycle is detected

The detection uses the community-total biomass, on a log scale and with
its mean removed, as a scalar signal, sampled after every time step. At
least 20 steps are needed before any cycle can be reported.

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
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md),
which works out the eigenvalues of the linearised dynamics instead of
watching a trajectory.

The reported `period` is a multiple of `dt`, so it is only resolved to
that accuracy; reduce `dt` if you need the period more precisely.

## What you get back may not be a steady state

The stopping criterion is a proxy. It says that two states `t_per` years
apart differ by less than `distance_tol` on whatever scale the criterion
is measured on; it does not say that the state reached is a fixed point.
There are four ways the returned object can fail to be one:

- the run converged on its own scale while the biomasses are still
  visibly drifting (`termination = "distance_tolerance"`);

- the run reached `t_max` without converging
  (`termination = "time_limit"`);

- the run settled on a limit cycle (`termination = "cycle_detected"`),
  in which case the state stored is one point on that cycle;

- the run stopped because a species was going extinct
  (`termination = "extinction"`).

So treat the result as a claim to be checked rather than as a guarantee:

    attr(params, "convergence")$attractor  # "fixed_point", "limit_cycle" or NA
    attr(params, "convergence")$residual   # largest biomass drift, in 1/year
    isSteady(params)                       # TRUE if within tolerance
    summary(params)                        # includes the biomass-drift verdict
    plot(getSteadyResidual(params))        # which species, and at which sizes

`attractor` is the field that answers the question: it is
`"fixed_point"` only where the measured biomass drift is within
`residual_tol`, so it cannot be satisfied by a distance function that
has merely gone quiet. `termination` says how the run ended and
`converged` whether the solver met its own criterion; neither is a claim
about the state. The last line says *where* the model is not steady,
which is the one to reach for when it is not: a model that is off steady
state is usually off in one species or one part of the size range, and
the plot names it. See
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md)
for why the verdict is phrased in terms of biomass drift rather than the
largest per-capita rate.

The messages this function prints say the same thing — a converged run
whose biomasses are still moving reports the drift and adds "Reduce the
tolerance on the distance function to converge further." — but they are
suppressed by `info_level = 0`, so in a script the `"convergence"`
attribute is the reliable check.

Finally, a genuine fixed point need not be a *stable* one. Use
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
to find out, and `solver = "newton"` to converge onto a fixed point that
the dynamics themselves would run away from.

## See also

[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md),
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md),
[`isSteady()`](https://sizespectrum.org/mizer/reference/isSteady.md),
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md),
[`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md),
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
