# Project values for first time step of Euler method

This is an internal function used by the user-facing
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
function. It is of potential interest only to mizer extension authors.

## Usage

``` r
project_n(
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
  flux_limiter = "none"
)

project_n_no_diffusion(
  params,
  r,
  n,
  dt,
  a,
  b,
  S,
  idx,
  w_min_idx_array_ref,
  no_sp,
  no_w
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

- flux_limiter:

  Name of the flux limiter used for a deferred high-order correction of
  the upwind advective flux, or `"none"` for plain first-order upwind.
  See
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

## Value

The updated abundance density matrix `n`.

## Details

The function calculates the abundance at the next time step using the
McKendrick-von Foerster equation: \$\$\frac{\partial N}{\partial t} +
\frac{\partial}{\partial w} \left( g N - \frac{1}{2}\frac{\partial(D
N)}{\partial w} \right) = -\mu N\$\$ which is solved using a
semi-implicit upwind finite volume scheme.

## See also

[`project`](https://sizespectrum.org/mizer/reference/project.md),
[`mizerRates`](https://sizespectrum.org/mizer/reference/mizerRates.md)
