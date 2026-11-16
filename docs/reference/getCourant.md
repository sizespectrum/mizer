# Growth Courant numbers on the consumer size grid

Measures how far growth can transport fish through a size bin during
`dt`. Each value is `dt / dw` times the larger positive growth velocity
at the bin's inflow and outflow faces. The face convention follows the
configured flux scheme: upwind uses `g[j - 1]` and `g[j]`, whereas the
second-order schemes use `g[j]` and `g[j + 1]`. Recruitment replaces
growth inflow at the egg boundary. Values outside the species' retained
size support are zero; empty but reachable bins are included.

## Usage

``` r
getCourant(
  params,
  dt = 1,
  n = initialN(params),
  n_pp = initialNResource(params),
  n_other = initialNOther(params),
  effort = params@initial_effort,
  t = 0,
  e_growth = NULL,
  ...
)
```

## Arguments

- params:

  A `MizerParams` object.

- dt:

  Positive time interval in years. Default 1 gives the Courant number
  per year, which scales linearly with the chosen interval.

- n:

  Consumer number densities. Defaults to the initial state.

- n_pp:

  Resource number densities. Defaults to the initial state.

- n_other:

  Other component abundances. Defaults to the initial state.

- effort:

  Fishing effort, as in
  [`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md).

- t:

  Time at which to calculate the rates, in years.

- e_growth:

  Optional precomputed growth matrix, with the same dimensions and, when
  present, species and size labels as `initialN(params)`. Must be
  finite. If omitted, growth is obtained from
  [`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md)
  at the supplied state, respecting custom rate functions.

- ...:

  Further arguments passed to the rate functions.

## Value

An `ArraySpeciesBySize` containing dimensionless growth Courant numbers.
With no positive growth every entry is zero.

## Details

A large value can indicate cell-scale oscillations with a high-order
flux. The target used by
[`project()`](https://sizespectrum.org/mizer/reference/project.md) is an
advisory growth-based measure, not a universal stability or positivity
bound. It does not measure mortality, diffusion or changes in rates
during the step. The exact bin widths are used, not the approximation
`w * log(w[2] / w[1])`.

## Examples

``` r
max(getCourant(NS_params, dt = 1 / 52))
#> [1] 12.70361
```
