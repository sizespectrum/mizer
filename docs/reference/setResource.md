# Set resource dynamics

Sets the intrinsic resource birth rate and the intrinsic resource
carrying capacity as well as the name of the function used to simulate
the resource dynamics. By default, the birth rate and the carrying
capacity are changed together in such a way that the resource
replenishes at the same rate at which it is consumed. So you should only
provide either the `resource_rate` or the `resource_capacity` (or
`resource_level`) because the other is determined by the requirement
that the resource replenishes at the same rate at which it is consumed.

## Usage

``` r
setResource(
  params,
  resource_rate = NULL,
  resource_capacity = NULL,
  resource_level = NULL,
  resource_dynamics = NULL,
  lambda = resource_params(params)[["lambda"]],
  n = resource_params(params)[["n"]],
  w_pp_cutoff = resource_params(params)[["w_pp_cutoff"]],
  balance = NULL,
  reset = FALSE,
  ...
)

resource_rate(params)

resource_rate(params, balance = NULL) <- value

resource_capacity(params)

resource_capacity(params, balance = NULL) <- value

resource_level(params)

resource_level(params, balance = NULL) <- value

resource_dynamics(params)

resource_dynamics(params, balance = NULL) <- value
```

## Arguments

- params:

  A MizerParams object

- resource_rate:

  Optional. A vector of per-capita resource birth rate for each size
  class or a single number giving the coefficient in the power-law for
  this rate, see "Setting resource dynamics" below. Must be strictly
  positive.

- resource_capacity:

  Optional. Vector of resource intrinsic carrying capacities or
  coefficient in the power-law for the capacity, see "Setting resource
  dynamics" below. The resource capacity must not be smaller than the
  resource abundance.

- resource_level:

  Optional. The ratio between the current resource number density and
  the resource capacity. Either a number used at all sizes or a vector
  specifying a value for each size. Must be greater than 0 and at most
  1, except at sizes where the resource is zero, where it can be `NaN`.
  This determines the resource capacity, so do not specify both this and
  `resource_capacity`.

- resource_dynamics:

  Optional. Name of the function that determines the resource dynamics
  by calculating the resource spectrum at the next time step from the
  current state.

- lambda:

  Used to set power-law exponent for resource capacity if the
  `resource_capacity` argument is given as a single number.

- n:

  Used to set power-law exponent for resource rate if the
  `resource_rate` argument is given as a single number.

- w_pp_cutoff:

  The upper cut off size of the resource spectrum power law used when
  `resource_capacity` is given as a single number. When changing
  `w_pp_cutoff` without providing `resource_capacity`, the cutoff can
  only be decreased. In that case, both the carrying capacity and the
  initial resource abundance will be cut off at the new value. To
  increase the cutoff, you must also provide the `resource_capacity` for
  the extended range.

- balance:

  By default, if possible, the resource parameters are set so that the
  resource replenishes at the same rate at which it is consumed. In this
  case you should only specify either the resource rate or the resource
  capacity (or resource level) because the other is then determined
  automatically. Set to FALSE if you do not want the balancing.

- reset:

  If set to TRUE, then the resource capacity and birth rate will be
  reset to the values calculated from the resource parameters, even if
  they were previously overwritten with custom values. If set to FALSE
  (default) then a recalculation from the resource parameters will take
  place only if no custom values have been set.

- ...:

  Unused

- value:

  The desired new value for the respective parameter.

## Value

`setResource`: A MizerParams object with updated resource parameters

A vector with the intrinsic resource birth rate for each size class.

A vector with the intrinsic resource capacity for each size class.

A vector with the ratio between the current resource number density and
the resource capacity for each size class.

The name of the function that determines the resource dynamics.

## Details

You would usually set the resource dynamics only after having finished
the calibration of the steady state. Then setting the resource dynamics
with this function will preserve that steady state, unless you
explicitly choose to set `balance = FALSE`. Your choice of the resource
dynamics only affects the dynamics around the steady state. The higher
the resource rate or the lower the resource capacity the less sensitive
the model will be to changes in the competition for resource.

If you provide the `resource_level` then that sets the
`resource_capacity` to the current resource number density divided by
the resource level. So in that case you should not specify
`resource_capacity` as well.

If you provide none of the arguments `resource_level`, `resource_rate`
or `resource_capacity`, and you do not change any of the resource
parameters, then the resource rate is kept at its previous value and,
when balancing, the capacity is recalculated from it. If instead you
change one of the resource parameters (`kappa`, `lambda`, `n` or
`w_pp_cutoff`) or set `reset = TRUE`, the rate and capacity are
recalculated from the resource parameters (and then balanced, unless
`balance = FALSE`).

## Setting resource dynamics

The `resource_dynamics` argument allows you to choose the resource
dynamics function. By default, mizer uses a semichemostat model to
describe the resource dynamics in each size class independently. This
semichemostat dynamics is implemented by the function
[`resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.md).
You can change that to use a logistic model implemented by
[`resource_logistic()`](https://sizespectrum.org/mizer/reference/resource_logistic.md)
or you can use
[`resource_constant()`](https://sizespectrum.org/mizer/reference/resource_constant.md)
which keeps the resource constant or you can write your own function.

Both the
[`resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.md)
and the
[`resource_logistic()`](https://sizespectrum.org/mizer/reference/resource_logistic.md)
dynamics are parametrised in terms of a size-dependent birth rate
\\r_R(w)\\ and a size-dependent capacity \\c_R\\. The help pages of
these functions give the details.

The `resource_rate` argument can be a vector (with the same length as
`w_full(params)`) specifying the intrinsic resource birth rate for each
size class. Alternatively it can be a single number that is used as the
coefficient in a power law: then the intrinsic birth rate \\r_R(w)\\ at
size \\w\\ is set to \$\$r_R(w) = r_R w^{n-1}.\$\$ The power-law
exponent \\n\\ is taken from the `n` argument.

The `resource_capacity` argument can be a vector specifying the
intrinsic resource carrying capacity for each size class. Alternatively
it can be a single number that is used as the coefficient in a truncated
power law: then the intrinsic carrying capacity \\c_R(w)\\ at size \\w\\
is set to \$\$c_R(w) = c_R\\ w^{-\lambda}\$\$ for all \\w\\ less than
`w_pp_cutoff` and zero for larger sizes. The power-law exponent
\\\lambda\\ is taken from the `lambda` argument.

The values for `lambda`, `n` and `w_pp_cutoff` are stored in a list in
the `resource_params` slot of the MizerParams object so that they can be
re-used automatically in the future. If you specify `resource_rate` or
`resource_capacity` as a single number, that coefficient is likewise
stored, as `r_pp` and `kappa` respectively. That list can be accessed
with
[`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md).

## See also

[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md)

## Examples

``` r
params <- NS_params
resource_dynamics(params)
#> [1] "resource_semichemostat"
resource_dynamics(params) <- "resource_constant"
```
