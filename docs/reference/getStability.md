# Analyse the dynamic stability of a mizer steady state

**\[experimental\]** Computes the eigenvalues of the linearised
one-step-ahead map at the steady state stored in `params@initial_n`.
These eigenvalues determine whether the steady state is dynamically
stable and, when a Hopf bifurcation is approached, the period of the
emergent limit cycle.

## Usage

``` r
getStability(
  params,
  reproduction = c("fixed", "dynamic"),
  effort = params@initial_effort,
  include_resource = FALSE,
  extinction_floor = 1e-06,
  h = 1e-04,
  dt = 1
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object whose `initial_n` holds the steady state to analyse. Typically
  the output of
  [`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md).

- reproduction:

  Whether the reproduction rate is held fixed (`"fixed"`, default) or
  run dynamically (`"dynamic"`) during the one-step evaluation. Must
  match the choice used when the steady state was computed.

- effort:

  The fishing effort to use. By default the initial effort stored in
  `params`.

- include_resource:

  If `FALSE` (default) the resource is treated as a quasi-static fast
  variable: for each perturbed fish abundance the resource is set to its
  analytic steady-state value conditioned on that fish abundance (valid
  only for semichemostat resource dynamics). If `TRUE`, both fish and
  resource cells are perturbed independently and the resource is evolved
  with the full resource dynamics function stored in
  `params@resource_dynamics`, giving the complete coupled Jacobian.

- extinction_floor:

  Relative abundance floor for the dynamic reproduction case. Default is
  `1e-6`.

- h:

  Relative step size for centred finite differences. Default `1e-4`. The
  result should not depend on this choice. If it does, the one-step map
  is not smooth at the state being analysed — see the section below.

- dt:

  The time step size to use for evaluating the numerical one-step map.
  Default `1`. The continuous eigenvalues are independent of this
  choice, but the discrete eigenvalues and `spectral_radius` returned
  will reflect the stability of mizer's numerical Euler method exactly
  for this step size.

## Value

A named list with the following components:

- `eigenvalues`:

  Complex vector of the valid continuous-time eigenvalues (\\\lambda_i =
  (1 - 1/\mu_i) / dt\\), sorted by decreasing real part. These describe
  the stability of the underlying continuous ODEs/PDEs.

- `discrete_eigenvalues`:

  Complex vector of the raw discrete eigenvalues \\\mu_i\\ of the
  numerical one-step map (evaluated at step size `dt`), sorted by
  decreasing real part of their continuous counterparts.

- `spectral_radius`:

  The spectral radius of the numerical one-step map evaluated at step
  size `dt`: \\\max_i\|\mu_i\|\\. A value less than 1 indicates that the
  numerical scheme is stable.

- `max_real_part`:

  The largest real part of the continuous eigenvalues: \\\max_i
  \text{Re}(\lambda_i)\\. Greater than 0 means unstable.

- `stable`:

  Logical: `TRUE` when `max_real_part < 0`.

- `dominant_period`:

  The period (in years) of the dominant continuous eigenvalue:
  `2*pi / abs(Im(lambda_1))`. `Inf` for a real positive dominant
  eigenvalue (monotone dynamics).

- `hopf_period`:

  Period (in years) of the complex continuous eigenvalue with the
  largest real part; `NULL` when no complex eigenvalue exists. This is
  the expected limit-cycle period near a Hopf bifurcation.

- `n_active`:

  Dimension of the Jacobian: number of active fish cells when
  `include_resource = FALSE`, or fish cells plus all resource cells when
  `include_resource = TRUE`.

- `leading_eigenvectors`:

  The eigenvectors of the two largest-modulus eigenvalues, reshaped back
  into the fish abundance space. When `include_resource = FALSE`: a
  complex array of shape `(n_species, n_sizes, 2)` with the same species
  and size dimnames as `params@initial_n`. When
  `include_resource = TRUE`: a list with `$fish` (the same array) and
  `$resource` (a complex matrix of shape `(n_w_full, 2)` for the
  resource component). Each eigenvector is normalised so that its
  maximum modulus equals 1. The real and imaginary parts of eigenvector
  1 span the two-dimensional oscillation plane of the dominant mode;
  [`Mod()`](https://rdrr.io/r/base/complex.html) gives the amplitude
  pattern across species and sizes.

## Details

### Mathematical background

The mizer time step applies a backward-Euler transport solve for the
fish: \$\$A(N^t, n\_{pp}^t)\\N^{t+1} = S(N^t, n\_{pp}^t),\$\$ and an
exact semi-chemostat update for the resource: \$\$n\_{pp}^{t+1} =
n\_{pp}^\* + (n\_{pp}^t - n\_{pp}^\*)\\e^{-\mu^t\\dt},\$\$ where
\\n\_{pp}^\* = r\_{pp}\\c\_{pp}/\mu^t\\ is the resource steady state
conditioned on the mortality \\\mu^t\\ due to consumers at time \\t\\.
Note that this function evaluates the Jacobian of this specific
first-order backward-Euler time step, regardless of which `method` you
might later pass to
[`project()`](https://sizespectrum.org/mizer/reference/project.md).
However, it fully respects any higher-order spatial scheme configured
via
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md).

The stability is determined by the Jacobian of the full one-step-ahead
map \\G : (N, n\_{pp}) \mapsto (N^{t+1}, n\_{pp}^{t+1})\\ at the fixed
point.

When `include_resource = FALSE` (the default), the resource is treated
as a fast variable that adjusts *instantaneously* to the consumer
abundance: for each perturbed \\N\\, \\n\_{pp}\\ is set to its
quasi-static equilibrium \\n\_{pp}^\*(N)\\. The resulting reduced
Jacobian \\L\_{\text{red}}\\ has dimension equal to the number of active
fish cells. This is equivalent to projecting the full dynamics onto the
slow manifold \\n\_{pp} = n\_{pp}^\*(N)\\.

When `include_resource = TRUE`, both fish and resource cells are
perturbed independently and the full coupled Jacobian
\\L\_{\text{full}}\\ is returned. Its eigenvalues include both the slow
fish modes and a cluster of fast resource-relaxation modes (with modulus
\\e^{-\mu\\dt} \ll 1\\). Comparing the dominant eigenvalues of the two
analyses shows how much the quasi-static approximation affects the
stability conclusion.

The discrete eigenvalues \\\mu_i\\ of the numerical Jacobian are mapped
back to their exact continuous-time equivalents \\\lambda_i = (1 -
1/\mu_i) / dt\\ to remove the artificial temporal numerical diffusion
introduced by the backward Euler solver. The steady state is **stable**
when all continuous-time eigenvalues satisfy \\\text{Re}(\lambda_i) \<
0\\ and **unstable** when at least one exceeds 0.

A **Hopf bifurcation** occurs when a complex-conjugate pair of
eigenvalues crosses the imaginary axis, giving a limit-cycle period
\$\$T = \frac{2\pi}{\|\text{Im}(\lambda)\|} \text{ years.}\$\$

Both branches use the same `project_n_loop()` C++ Thomas solver as the
regular dynamics, evaluating the transport coefficients with the exact
spatial scheme configured in `params` (e.g., first-order upwind or a
second-order limiter). The Jacobian is computed numerically using a
multiplicative (relative) finite-difference step \\h \cdot N^\*\\. Where
a cell sits at exactly zero and so has no scale of its own, the step is
floored at the local scale of the spectrum, interpolated from the
nonzero neighbours, so that the cell still gets a resolved derivative
rather than a column of rounding error.

Every state at which the rate functions are evaluated satisfies \\N \ge
0\\: where a centred step would push a cell negative — which can only
happen for a cell at (or below) the floor described above — the column
is differenced forwards from the unperturbed state instead. At the
boundary of the physical cone the one-sided derivative is the
appropriate object anyway, since the dynamics never visit the states a
centred step would sample. A rate function registered with
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
therefore never has to be defined at negative abundances. Such columns
are first order in `h` rather than second, so they respond slightly more
to a change of `h` than the rest.

## Requires a smooth one-step map

The finite-difference Jacobian is only meaningful if the one-step map is
differentiable at \\N^\*\\. A custom rate function registered with
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
that jumps as a function of the abundances breaks this in two ways. If
the state sits on the switching threshold, some perturbations straddle
it and pick up the jump, and the reported spectral radius then varies
wildly with `h`. If the state is near but not on the threshold, no
perturbation crosses it, and the function silently returns the stability
of the single branch the state happens to lie on — which can read as
`stable` for a model whose simulations never settle.

Re-running with a different `h` is the cheapest check: if the answer
moves, do not trust it. See [Discontinuous rate
functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html).

## See also

[`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md)
