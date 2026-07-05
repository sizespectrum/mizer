# Project values with a predictor-corrector method

This is an experimental second-order time stepping variant of
[`project_n()`](https://sizespectrum.org/mizer/reference/project_n.md).
It first predicts the new consumer densities with
[`project_n()`](https://sizespectrum.org/mizer/reference/project_n.md),
optionally recalculates rates from that prediction, and then applies a
Crank-Nicolson corrector using midpoint rates.

## Usage

``` r
project_n_2(
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

  Optional midpoint rates. If supplied, these are used directly in the
  Crank-Nicolson corrector.

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

If the rate recalculation arguments are not supplied, the corrector uses
the supplied rates as fixed rates. In that case the corrector is second
order only for the frozen-rate transport problem, not for the full
nonlinear mizer dynamics.

## See also

[`project_n`](https://sizespectrum.org/mizer/reference/project_n.md)
