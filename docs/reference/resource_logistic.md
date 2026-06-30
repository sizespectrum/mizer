# Project resource using logistic model

If you set your resource dynamics to use this function then the time
evolution of the resource spectrum is described by a logistic equation
\$\$\frac{\partial N_R(w,t)}{\partial t} = r_R(w) N_R(w)\Big\[ 1 -
\frac{N_R(w,t)}{c_R (w)} \Big\] - \mu_R(w, t) N_R(w,t)\$\$

## Usage

``` r
resource_logistic(
  params,
  n,
  n_pp,
  n_other,
  rates,
  t,
  dt,
  resource_rate,
  resource_capacity,
  ...
)

balance_resource_logistic(params, resource_rate, resource_capacity)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.md)
  object

- n:

  A matrix of species abundances (species x size)

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list with the abundances of other components

- rates:

  A list of rates as returned by
  [`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md)

- t:

  The current time

- dt:

  Time step

- resource_rate:

  Resource replenishment rate

- resource_capacity:

  Resource carrying capacity

- ...:

  Unused

## Value

Vector containing the resource number density in each size class at the
next timestep

## Details

Here \\r_R(w)\\ is the resource regeneration rate and \\c_R(w)\\ is the
carrying capacity in the absence of predation. These parameters are
changed with
[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md).
The mortality \\\mu_R(w, t)\\ is due to predation by consumers and is
calculate with
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md).

This function uses the analytic solution of the above equation to
calculate the resource abundance at time `t + dt` from all abundances
and rates at time `t`, keeping the mortality fixed during the timestep.

To set your model to use logistic dynamics for the resource you do

    params <- setResource(params,
                          resource_dynamics = "resource_logistic",
                          resource_level = 0.5)

where you should replace `params` with the name of the variable holding
your MizerParams object. You can of course choose any value between 0
and 1 for the resource level.

The `balance_resource_logistic()` function is called by
[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md)
to determine the values of the resource parameters that are needed to
make the replenishment rate at each size equal the consumption rate at
that size, as calculated by
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md).
It should be called with exactly one of `resource_rate` or
`resource_capacity` and returns a named list with values for both. If
`resource_rate` is supplied it must be at least as large as the current
mortality at each size. If `resource_capacity` is supplied it must be
not be less than the current resource abundance. Where it equals the
current resource abundance and there is positive consumption, it is
nudged upwards slightly to avoid division by zero.

## See also

[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md)

Other resource dynamics functions:
[`resource_constant()`](https://sizespectrum.org/mizer/reference/resource_constant.md),
[`resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.md)
