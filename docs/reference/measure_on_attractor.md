# Measure a quantity on the attractor a projection settled on

Measure a quantity on the attractor a projection settled on

## Usage

``` r
measure_on_attractor(
  settled,
  value_func,
  conv,
  dt,
  t_sample,
  sample_all,
  method,
  default_name = "Value"
)
```

## Arguments

- settled:

  The MizerParams returned by `project_until_settled()`.

- value_func:

  The function measuring the quantity.

- conv:

  The `"convergence"` attribute of `settled`.

- dt:

  The time step.

- t_sample:

  The averaging window to use when nothing settled.

- sample_all:

  Whether to sample even at a fixed point.

- method:

  The numerical method.

- default_name:

  The series name to use when `value_func` supplies none.

## Value

A list with the mean, minimum and maximum over the attractor, the names
of the series and the metadata read off `value_func`'s result.
