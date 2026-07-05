# Project values with the TR-BDF2 method

This is an L-stable, second-order time stepping variant of
[`project_n()`](https://sizespectrum.org/mizer/reference/project_n.md).
It takes one TR-BDF2 step, consisting of a trapezoidal (Crank-Nicolson)
stage over the first part of the time step followed by a second-order
backward differentiation (BDF2) stage over the remainder.

## Usage

``` r
project_n_tr_bdf2(
  params,
  r,
  n,
  dt,
  a,
  b,
  c,
  S,
  idx,
  w_min_idx_array_ref,
  no_sp,
  no_w,
  rates_fns = NULL,
  n_pp = NULL,
  n_other = NULL,
  t = 0,
  effort = NULL,
  r_hat = NULL,
  r_mid = NULL,
  n_hat = NULL,
  flux_limiter = "none",
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object.

- r:

  A list of rates as returned by
  [`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md).

- n:

  An array (species x size) with the number density at the current time
  step.

- dt:

  Time step.

- a:

  A matrix (species x size) used in the solver (transport term).

- b:

  A matrix (species x size) used in the solver (diagonal term).

- c:

  A matrix (species x size) used in the solver (transport term).

- S:

  A matrix (species x size) used in the solver (source term).

- idx:

  Index vector for size bins (excluding the first one).

- w_min_idx_array_ref:

  Index vector for the start of the size spectrum for each species.

- no_sp:

  Number of species.

- no_w:

  Number of size bins.

- rates_fns:

  Optional named list of rate functions, as used by
  [`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md).
  If supplied together with `n_pp`, `n_other` and `effort`, provisional
  end-of-step rates are calculated from the predicted densities.

- n_pp:

  Resource abundance used when recalculating provisional rates.

- n_other:

  Other ecosystem components used when recalculating provisional rates.

- t:

  Current time.

- effort:

  Fishing effort used when recalculating provisional rates.

- r_hat:

  Optional provisional end-of-step rates. If supplied, these are used
  instead of recalculating them.

- r_mid:

  Optional midpoint rates. If supplied, these are used directly to build
  the TR-BDF2 operator.

- n_hat:

  Optional provisional end-of-step densities. When supplied with a flux
  limiter, the limiter is frozen at the midpoint field
  `(n + n_hat) / 2`; otherwise it falls back to the start-of-step `n`.

- flux_limiter:

  Name of the flux limiter used for a deferred high-order correction of
  the upwind advective flux, or `"none"` for plain first-order upwind.
  See
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- ...:

  Further arguments passed to the rate functions.

## Value

The updated abundance density matrix `n`.

## Details

The nonlinear rates are handled exactly as in
[`project_n_2()`](https://sizespectrum.org/mizer/reference/project_n_2.md):
a provisional Euler predictor gives end-of-step rates, which are
averaged with the start-of-step rates to obtain second-order-accurate
midpoint rates `r_mid`. Both TR-BDF2 stages then use this single frozen
operator.

With the standard parameter \\\gamma = 2 - \sqrt 2\\ both stages share
the same implicit coefficient \\\alpha\\\Delta t\\ with \\\alpha =
\gamma/2 = 1 - 1/\sqrt 2\\, so the operator \\I - \alpha\\\Delta t\\L\\
is assembled once with `get_transport_coefs()` and each stage is a
single tridiagonal solve with `project_n_loop()`, exactly as in
[`project_n()`](https://sizespectrum.org/mizer/reference/project_n.md)
and
[`project_n_2()`](https://sizespectrum.org/mizer/reference/project_n_2.md).
Unlike the Crank-Nicolson corrector in
[`project_n_2()`](https://sizespectrum.org/mizer/reference/project_n_2.md),
TR-BDF2 is L-stable and therefore damps the stiff modes that cause
Crank-Nicolson to oscillate at large time steps.

If the rate recalculation arguments are not supplied, the step uses the
supplied rates as fixed rates. In that case the method is second order
only for the frozen-rate transport problem, not for the full nonlinear
mizer dynamics, but it remains L-stable.

## See also

[`project_n`](https://sizespectrum.org/mizer/reference/project_n.md),
[`project_n_2`](https://sizespectrum.org/mizer/reference/project_n_2.md)
