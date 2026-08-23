# Analyse the stability of mizer's numerical time step

**\[experimental\]** Computes the eigenvalues \\\mu_i\\ of the
linearised one-step-ahead map at the steady state stored in
`params@initial_n`, for a given step size `dt`. These describe how
mizer's numerical scheme, rather than the model, behaves near the steady
state: the map does not amplify perturbations when the spectral radius
\\\max_i\|\mu_i\|\\ is less than 1.

## Usage

``` r
getDiscreteStability(params, effort = params@initial_effort, h = 1e-04, dt = 1)
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

- dt:

  The time step size of the one-step map. Default `1`.

## Value

A named list with the following components:

- `discrete_eigenvalues`:

  Complex vector of the eigenvalues \\\mu_i\\ of the one-step map,
  sorted by decreasing modulus.

- `spectral_radius`:

  \\\max_i\|\mu_i\|\\. Less than 1 means the numerical scheme is stable
  at this `dt`.

- `stable`:

  Logical: `TRUE` when `spectral_radius < 1`.

- `dt`:

  The step size the map was evaluated at.

- `n_active`:

  Dimension of the Jacobian.

- `leading_eigenvectors`:

  The eigenvectors of the two largest-modulus eigenvalues, in the same
  shape as for
  [`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md).

- `params`:

  The validated `params` object the analysis was made at.

## Details

This is the numerical counterpart of
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md),
which analyses the model itself and involves no time step at all. Use
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
to ask whether the steady state of the *model* is stable, and this
function to ask what mizer's solver does at a particular `dt`. The two
can disagree, and that disagreement is the point: the implicit transport
solve damps oscillations artificially, so a physically unstable steady
state can have a spectral radius below 1 at a large `dt`, and the
simulation then sits at a state the model does not actually hold.

### The map that is linearised

One step is what
[`project()`](https://sizespectrum.org/mizer/reference/project.md) takes
with `method = "euler"`: the rates are evaluated at the state at the
start of the step, and the resulting transport problem is solved
implicitly, \$\$A(N^t, n\_{pp}^t)\\N^{t+1} = S(N^t, n\_{pp}^t),\$\$ with
the same `project_n_loop()` C++ Thomas solver and the same spatial
scheme
([`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md))
as the regular dynamics.

Because the rates are evaluated at the *input* state, the step is not
fully implicit, and the discrete eigenvalues therefore cannot be
converted into continuous-time eigenvalues by any exact algebraic
relation. That conversion is what
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
avoids by differentiating the rates of change themselves.

The resource is advanced by the model's own `resource_dynamics`
function, the one
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
calls. Nothing is substituted for it: the map that is differentiated
here reproduces a single `project(method = "euler")` step exactly, which
is what makes the spectral radius a statement about mizer's solver
rather than about a nearby scheme.

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

[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md),
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
