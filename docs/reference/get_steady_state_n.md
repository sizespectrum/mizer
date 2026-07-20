# Calculate steady state abundance

This function calculates the steady state abundance by solving the
transport equation with given growth and mortality rates. It sets up a
tri-diagonal system and solves it.

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
  relax = 0.3
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

  Maximum number of Picard iterations used when a flux limiter is
  active.

- tol:

  Relative convergence tolerance for the Picard iteration.

- relax:

  Under-relaxation factor in (0, 1\] for the Picard iteration when a
  flux limiter is active.

## Value

A matrix with the steady state abundance

## Details

The spatial discretisation of the advective flux is read from the `flux`
entry of the `second_order_w` slot of `params`. With a second-order flux
scheme active the steady state must match the one that
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
converges to. Because the limiter depends on the solution, the steady
state is then found by an under-relaxed Picard iteration: the limiter is
frozen at the current iterate, the resulting tridiagonal system is
solved, and the iterate is updated towards that solution, repeating
until it converges. (At `dt = 1` the limited operator is not diagonally
dominant, so the plain fixed-point map only stalls; under-relaxation
makes it converge.)

The returned abundance is held at zero above each species' `w_max`, the
same upper boundary condition that
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
imposes (via `zero_above_support()` in
[`project_n()`](https://sizespectrum.org/mizer/reference/project_n.md)).
Without this the bottom-up solve would carry density above `w_max`
whenever growth is still positive there or diffusion pushes density past
it.
