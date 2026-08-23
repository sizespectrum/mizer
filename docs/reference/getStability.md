# Analyse the dynamic stability of a mizer steady state

**\[experimental\]** Computes the eigenvalues of the linearised dynamics
at the steady state stored in `params@initial_n`. These eigenvalues
determine whether the steady state is dynamically stable and, where the
spectrum contains a complex pair, the period at which the model
oscillates.

## Usage

``` r
getStability(params, effort = params@initial_effort, h = 1e-04)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object whose `initial_n` holds the steady state to analyse. Typically
  the output of
  [`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md).

- effort:

  The fishing effort to use. By default the initial effort stored in
  `params`.

- h:

  Relative step size for centred finite differences. Default `1e-4`. The
  result should not depend on this choice. If it does, the dynamics are
  not smooth at the state being analysed — see the section below.

## Value

A named list with the following components:

- `eigenvalues`:

  Complex vector of the continuous-time eigenvalues \\\lambda_i\\,
  sorted by decreasing real part.

- `max_real_part`:

  The largest real part of the eigenvalues: \\\max_i
  \text{Re}(\lambda_i)\\. Greater than 0 means unstable.

- `stable`:

  Logical: `TRUE` when `max_real_part < 0`.

- `dominant_period`:

  The period (in years) of the dominant eigenvalue:
  `2*pi / abs(Im(lambda_1))`. `Inf` for a real dominant eigenvalue
  (monotone dynamics).

- `oscillation_period`:

  Period (in years) of the oscillatory mode below, \\2\pi/\|\omega\|\\;
  `NULL` when no complex eigenvalue exists. It is the period at which
  the model rings, and only at a Hopf bifurcation — where the real part
  is zero — the period of a limit cycle.

- `leading_oscillatory_eigenvalue`:

  The complex eigenvalue with the largest real part, or `NULL` when
  there is none. Its real part is the rate at which that oscillation
  grows, so a strongly negative one means the ringing is a transient,
  not a cycle the model settles onto.

- `leading_oscillatory_eigenvector`:

  Its eigenvector, as a list with `$fish`, a complex
  `(n_species, n_sizes)` matrix, and `$resource`, a complex vector of
  length `n_w_full`. This is the mode
  [`getOscillationModeSim()`](https://sizespectrum.org/mizer/reference/getOscillationModeSim.md)
  draws, and it is not in general one of `leading_eigenvectors`: the
  dominant mode of the system can be real while the dominant
  *oscillatory* mode is well down the spectrum.

- `n_active`:

  Dimension of the Jacobian: the number of active fish cells plus all
  resource cells.

- `leading_eigenvectors`:

  The eigenvectors of the two eigenvalues with the largest real part,
  reshaped back into the state space: a list with `$fish`, a complex
  array of shape `(n_species, n_sizes, 2)` with the same species and
  size dimnames as `params@initial_n`, and `$resource`, a complex matrix
  of shape `(n_w_full, 2)`. Each eigenvector is normalised by a single
  scalar covering both blocks, so that the relative amplitude and phase
  between fish and resource are those of the mode. The scalar is chosen
  so that the largest perturbation *relative to the steady state*,
  \\\|v_i\| / x^\*\_i\\, is 1 somewhere in the state: an absolute
  normalisation would be set entirely by the resource, whose densities
  dwarf the fish abundances. `Mod(fish[, , 1]) / initialN(params)` is
  therefore the relative amplitude pattern, peaking at 1 in whichever
  cell swings hardest. The real and imaginary parts of eigenvector 1
  span the two-dimensional oscillation plane of the dominant mode.

- `params`:

  The validated `params` object the analysis was made at.

## Details

### Mathematical background

Mizer discretises the size axis but not time: on the size grid the model
is a system of ordinary differential equations \$\$\frac{dN}{dt} = F(N,
n\_{pp}),\$\$ where \\F\\ collects the divergence of the growth flux,
the mortality sink and the reproductive influx at the egg size,
assembled with the spatial scheme configured via
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md).
`getStability()` differentiates \\F\\ directly, by centred finite
differences in each state variable, and returns the eigenvalues
\\\lambda_i\\ of the resulting Jacobian \\J = \partial F/\partial N\\.
The steady state is **stable** when all of them satisfy
\\\text{Re}(\lambda_i) \< 0\\ and **unstable** when at least one exceeds
0.

No time step enters this calculation. The eigenvalues are a property of
the model, not of any solver: they describe the continuous-time dynamics
of the semi-discretised model, and are what a simulation with a small
enough time step converges to. The stability of the numerical step
itself is a separate question, answered by
[`getDiscreteStability()`](https://sizespectrum.org/mizer/reference/getDiscreteStability.md).

A complex-conjugate pair \\\lambda = \sigma \pm i\omega\\ is an
oscillatory mode: a perturbation along it rings with period \$\$T =
\frac{2\pi}{\|\omega\|} \text{ years,}\$\$ growing or decaying as
\\e^{\sigma t}\\. The pair with the largest \\\sigma\\ is returned as
`leading_oscillatory_eigenvalue`, with its period and eigenvector.

That is a statement about the mode, not about a bifurcation. A **Hopf
bifurcation** is the event of such a pair *crossing* the imaginary axis,
and a single spectrum cannot show a crossing: the leading oscillatory
mode of a comfortably stable model can sit far to the left, ringing only
as a transient on the way back to the fixed point. Establishing a Hopf
bifurcation means watching \\\sigma\\ pass through zero as a parameter
is varied, which is what
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
is for. Only then is \\T\\ the period of an emerging limit cycle;
otherwise it is the period of a damped (or growing) oscillation.

### What is in the Jacobian

The resource is a state variable of the system like any other: fish and
resource cells are perturbed independently, giving the full coupled
Jacobian. Its eigenvalues include both the slow fish modes and a cluster
of fast resource-relaxation modes, at \\\lambda \approx -(r\_{pp} +
\mu_R)\\. Any resource dynamics function is supported: the semichemostat
derivative is written down analytically, and anything else is
differenced over a short step.

Components registered with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
are *not* state variables here. They are held at their stored values
while the fish and the resource are perturbed, so the spectrum is that
of the consumer-resource subsystem with the components frozen. This is
exact when a component is a fixed input, and a good approximation when
it is much faster or much slower than the fish, but it is not the full
model, and mizer says so with a warning when it meets one. Giving
extension components an explicit residual and Jacobian is the work that
would lift this restriction.

Reproduction is a state-dependent rate like any other: the reproduction
function stored in `params@rates_funcs$RDD` is evaluated at each
perturbed state, so the feedback from the spectra back onto the influx
of eggs is part of the Jacobian, exactly as it is part of
[`project()`](https://sizespectrum.org/mizer/reference/project.md).
There is no option to pin the reproduction rate at its value at the
fixed point. A model in which reproduction really is constant expresses
that as a model: with `rates_funcs$RDD = "constantRDD"` the derivative
of the reproduction rate is zero and the pinned Jacobian is what the
analysis returns.

This is why the stability of a steady state depends on the reproduction
parameters even though the steady state itself does not.
[`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
moves along a family of `erepro`/`R_max` pairs that all leave the same
fixed point, but they do not all leave the same dynamics: at a
[`reproduction_level()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
near 1 the reproduction rate barely responds to the energy invested in
it, approaching the constant-reproduction case, while at a level near 0
it follows that energy proportionally. The two ends can differ in their
verdict, so the analysis has to read the model rather than take an
argument.

### Numerical details

The Jacobian is computed numerically using a multiplicative (relative)
finite-difference step \\h \cdot N^\*\\. Where a cell sits at exactly
zero and so has no scale of its own, the step is floored at the local
scale of the spectrum, interpolated from the nonzero neighbours, so that
the cell still gets a resolved derivative rather than a column of
rounding error.

Every state at which the rates are evaluated satisfies \\N \ge 0\\:
where a centred step would push a cell negative — which can only happen
for a cell at (or below) the floor described above — the column is
differenced forwards from the unperturbed state instead. At the boundary
of the physical cone the one-sided derivative is the appropriate object
anyway, since the dynamics never visit the states a centred step would
sample. A rate function registered with
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
therefore never has to be defined at negative abundances. Such columns
are first order in `h` rather than second, so they respond slightly more
to a change of `h` than the rest.

## Requires smooth dynamics

The finite-difference Jacobian is only meaningful if the rates of change
are differentiable at \\N^\*\\. A custom rate function registered with
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
that jumps as a function of the abundances breaks this in two ways. If
the state sits on the switching threshold, some perturbations straddle
it and pick up the jump, and the reported eigenvalues then vary wildly
with `h`. If the state is near but not on the threshold, no perturbation
crosses it, and the function silently returns the stability of the
single branch the state happens to lie on — which can read as `stable`
for a model whose simulations never settle.

Re-running with a different `h` is the cheapest check: if the answer
moves, do not trust it. See [Discontinuous rate
functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html).

## See also

[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md),
[`getDiscreteStability()`](https://sizespectrum.org/mizer/reference/getDiscreteStability.md),
[`getOscillationModeSim()`](https://sizespectrum.org/mizer/reference/getOscillationModeSim.md)
