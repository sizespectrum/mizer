# Calculate steady state abundance

Solves the steady-state transport equation with growth, mortality and
diffusion held fixed and the supplied abundance in each egg size class.

## Usage

``` r
get_steady_state_n(
  params,
  g,
  mu,
  D,
  N0,
  max_iterations = 500,
  tol = 1e-10,
  relax = 0.5
)
```

## Arguments

- params:

  A MizerParams object

- g:

  A matrix of growth rates (species x size)

- mu:

  A matrix of mortality rates (species x size)

- D:

  A matrix of diffusion rates (species x size)

- N0:

  A vector with the abundance at the smallest size for each species

- max_iterations:

  Maximum number of defect-correction updates for a van Leer flux
  limiter.

- tol:

  Relative tolerance for both the correction and the transport residual,
  checked separately for every species.

- relax:

  Under-relaxation factor in (0, 1\] for defect correction.

## Value

A matrix with the steady-state abundance. Stops with an error naming the
unconverged species if the iteration limit is reached.

## Details

The flux scheme is read from
[`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md).
The upwind and centred systems are linear and are solved directly. With
a van Leer limiter, defect correction solves the nonlinear equation
using an upwind preconditioner with the same face placement and
diffusion as the second-order scheme. Only the preconditioner has its
reconstruction weights set to zero: the residual retains the model's
chosen flux limiter.

Each present species is solved with unit egg abundance and then rescaled
to `N0`, because the limiter is unchanged by multiplying a species'
spectrum by a positive constant. This prevents species with large egg
abundances from hiding non-convergence in other species. Convergence
requires both a small correction and a small residual of the original
transport equation, so clipping negative iterates cannot turn a stalled
iteration into a solution.

Abundance is held at zero below the egg size and above the support top
given by `support_top_idx()`, matching the boundary conditions of
[`project()`](https://sizespectrum.org/mizer/reference/project.md).
