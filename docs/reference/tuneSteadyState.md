# Tune a model so that the state it is in becomes a steady state

**\[experimental\]**

Solves for a model's steady state using declared tuning strategies. The
built-in defaults hold final recruitment and the resource abundance
fixed, then adjust reproduction parameters and resource capacity.
Extensions declare which states their strategies hold or evolve, what
parameters may change, and how the held targets are balanced. The
restored candidate is checked against every effective equation rate
before its final steady assessment.

It is not the function for asking where a model with given parameters
settles. Because the reproduction rate is held at its current value, it
cannot reach a steady state in which a species' reproduction rate is
very different, such as an alternative stable state or a state that
exists because predation changes the reproduction of the prey. Because
it re-derives `erepro`/`R_max` and `cc_pp`, two states found with it do
not share parameters. Use
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
for those questions.

## Usage

``` r
tuneSteadyState(
  params,
  solver = c("project", "newton"),
  effort = params@initial_effort,
  preserve = c("reproduction_level", "erepro", "R_max"),
  info_level = default_info_level(),
  scope = c("whole", "conditional"),
  fixed = character(),
  criterion = c("state", "biomass"),
  tolerances = NULL,
  atol = 1e-12,
  strategies = NULL,
  joint = NULL,
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- solver:

  The solver to use: `"project"` to run the dynamics until they settle,
  `"newton"` to solve the steady-state equation directly. See *Choosing
  a solver*.

- effort:

  The fishing effort to use throughout. By default the initial effort
  stored in `params`.

- preserve:

  **\[experimental\]** An explicit reproduction-preservation choice,
  such as `"erepro"`, `"R_max"` or `"reproduction_level"` for
  Beverton-Holt reproduction. Omitting this argument selects the
  provider's declared default strategy. Other laws support their own
  controls; unsupported choices error. Do not combine this selector with
  an explicit reproduction entry in `strategies`.

- info_level:

  Controls the amount of information messages that are shown. Higher
  levels lead to more messages, `info_level = 0` gives silence. The
  default is taken from the `mizer_info_level` option, see
  [`default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.md).

- scope:

  `"whole"` or `"conditional"`.

- fixed:

  Owners held fixed for a conditional assessment. Use `"consumers"`,
  `"resource"`, or `"component:name"` for a component named `name`. This
  must be empty for whole-model assessments and nonempty for conditional
  assessments.

- criterion:

  `"state"` for the absolute state residual, or `"biomass"` for net
  biomass drift in the built-in checks.

- tolerances:

  Optional named numeric vector overriding individual tolerances. Names
  are `owner/measure`, as shown in the `measures` table.

- atol:

  Absolute tolerance for built-in checks, in grams per year.

- strategies:

  Optional named character vector selecting a strategy for each
  specified owner, such as `c(resource = "evolve")` or
  `c("component:pool" = "hold")`. Omitted owners use their declared
  defaults.

- joint:

  Names of joint strategies registered with
  [`setTuningStrategy()`](https://sizespectrum.org/mizer/reference/setTuningStrategy.md).
  `NULL` selects declared joint defaults;
  [`character()`](https://rdrr.io/r/base/character.html) selects none. A
  joint strategy replaces the local defaults of its owners. Conflicting
  explicit local selections or overlapping joint ownership error.

- ...:

  Arguments for the chosen solver.

  With `solver = "project"`: `t_max`, `t_check`, `dt`, `distance_tol`,
  `residual_tol`, `amplitude_tol`, `amp_rel_tol`,
  `extinction_threshold`, `progress_bar`, `method`, `subcycle`,
  `max_courant` and `max_substeps`, all as described in
  [`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md).
  There is no `t_save`, because no trajectory is returned. Note that
  `distance_tol` here defaults to `0.1 * dt` and measures the largest
  relative change in egg production, because the distance function is
  [`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md).
  `residual_tol` is judged on the model as the search sees it, with
  reproduction, the resource and the other components pinned; the
  residual reported in the result is measured again on the model that is
  actually returned.

  With `solver = "newton"`: `solver_tol` (default `1e-6`), a tolerance
  on the per-capita rate of change passed to
  [`nleqslv::nleqslv()`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html).
  It was called `residual_tol` before mizer 3.3, a name now used by the
  assessment criterion above; `maxit` (default `200`); `jacobian`,
  either `"update"` (default, the Jacobian is computed once and then
  updated cheaply each iteration — `nleqslv`'s `"Broyden"`) or
  `"recompute"` (a numerical Jacobian at every iteration — `nleqslv`'s
  `"Newton"`); `global`, the globalisation strategy (default `"dbldog"`,
  a robust double-dogleg trust region); and `verbose` to trace the
  iterations.

## Value

A `MizerParams` object with the initial state replaced by the steady
state found, and with `erepro`/`R_max` and `cc_pp` adjusted as described
above. It carries a `"convergence"` attribute describing the solution
found; see
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md).
Check it: convergence is not guaranteed.

## Details

Each provider supplies named strategies and a default. Use `strategies`
to choose alternatives, and `joint` for declared strategies owning
coupled processes. The selected strategies must cover the whole tuning
operation without overlapping parameter ownership. Missing support,
unknown stationary behaviour and unsupported solver combinations stop
before the search.

Recruitment targets come from the complete effective rate pipeline,
including named replacements and extension wrappers. A held target
replaces the final recruitment entering the equations, so wrappers are
not applied a second time to a total that already includes their
contribution.

Balancers run once after the constrained search. They must preserve the
candidate state and all effective target equation rates, including
growth, diffusion, total mortality and final recruitment. Mizer
independently checks the released model; it does not repair changed
rates by solving the consumers again or by choosing an undeclared
numerical balancing method.

The `"convergence"` attribute separates `search_assessment`, measured
with the selected state and rate constraints, from `assessment`,
measured after restoring and balancing the model. Whole-model assessment
is the default. Select `scope = "conditional"` with explicit `fixed`
owners for a subsystem; those inputs must also be held by the selected
tuning strategies.

Holding those inputs fixed is what makes the search reliable, but it
does not make the result certain: the state that is stored is only as
close to a fixed point as the solver's tolerance allowed, and with
`solver = "project"` the run may instead have stopped on a limit cycle
or on a species going extinct. Check the result rather than assuming it;
see the section below.

## Choosing a solver

`solver = "project"` (the default) runs the dynamics until they settle,
via
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md),
using
[`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md)
as its distance function. It needs no extra packages. Every required
dynamic state must provide convergence support.

`solver = "newton"` solves the steady-state equation directly with a
Newton-type root finder from the `nleqslv` package. It converges even
when the steady state is dynamically **unstable**, where the
time-stepping solver cannot, and it discovers the support of the steady
state automatically. It starts from the spectra in `initialN(params)`,
so a reasonable initial guess still matters — for example the spectra
from a nearby stable parameterisation, or the (diverging) output of
`solver = "project"`.

Because the resource is held fixed either way, `solver = "newton"` here
does not need the resource to be a semichemostat, unlike in
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
where the resource is one of the unknowns. The Newton solver does not
solve extension component states; dynamic components require explicit
conditional exclusions with `fixed`, or use the projection solver.

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

Even a genuine fixed point can be one the model cannot sustain. `erepro`
is the proportion of the energy invested in reproduction that ends up as
egg biomass, so if a re-derived `erepro` is greater than 1, the
reproduction rate held during the search is inconsistent with the state
that was found. Usually this means that the state you asked for does not
exist with births held at their current values; try
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
instead. mizer warns when this happens, but with the default `preserve`
the warning reads the same for an `erepro` of 1.01 as for one of 1e34,
and `info_level = 0` suppresses it, so look at the values:

    species_params(params)$erepro          # greater than 1 is unrealistic

Finally, a genuine fixed point need not be a *stable* one. Use
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
to find out, and `solver = "newton"` to converge onto a fixed point that
the dynamics themselves would run away from.

## See also

[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md),
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md),
[`steadySingleSpecies()`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md),
[`isSteady()`](https://sizespectrum.org/mizer/reference/isSteady.md),
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)

## Examples

``` r
# \donttest{
params <- newTraitParams()
species_params(params)$gamma[5] <- 3000
params <- tuneSteadyState(params)
#> Reached the convergence tolerance after 21 years. The biomasses change at up to 0.0026 per year. This is a conditional assessment with resource, consumers/recruitment held fixed.
plotSpectra(params)

# }
```
