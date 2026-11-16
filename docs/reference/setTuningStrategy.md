# Register a joint tuning strategy

A joint strategy owns all listed responsibilities. Selecting it
suppresses their local defaults; overlapping owners or adjustable
parameter paths are errors. Its equation-rate callback returns one named
rate list per owner. The declaration is bound to the current effective
dispatch, including named rate replacements and extension methods.

## Usage

``` r
setTuningStrategy(
  params,
  name,
  owners,
  definition,
  equation_rates_fun,
  stationary,
  default = FALSE
)
```

## Arguments

- params:

  A MizerParams object.

- name:

  Unique name for the joint strategy.

- owners:

  Canonical tuning owners: consumers, resource, reproduction, or
  component:name.

- definition:

  A tuning strategy specification.

- equation_rates_fun:

  Callback returning effective equation quantities.

- stationary:

  Explicit declaration that this strategy supports stationary tuning of
  its processes. A prescribed state alone does not establish this.

- default:

  Whether tuning selects this joint strategy by default.

## Value

The updated model.
