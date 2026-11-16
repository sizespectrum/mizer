# Assess the steady state of the whole model or a conditional subsystem

**\[experimental\]** Returns a structured assessment of every dynamical
owner in the requested scope. The default includes the consumers,
resource and all components. Missing convergence support is reported
explicitly; it never counts as a successful check. A conditional
assessment requires an explicit `fixed` list and describes the subsystem
with those owners held at their current values. Owners declared
prescribed are recorded separately. Structural prescribed-coordinate
masks exclude only those cells from the assessment; the remaining cells
of a partially prescribed owner still require coverage.

## Usage

``` r
getSteadyAssessment(
  params,
  tol = 0.05,
  effort = params@initial_effort,
  scope = c("whole", "conditional"),
  fixed = character(),
  criterion = c("state", "biomass"),
  tolerances = NULL,
  atol = 1e-12,
  t = 0,
  dt = 1e-04,
  ...
)
```

## Arguments

- params:

  A MizerParams object.

- tol:

  Relative tolerance for built-in checks, in 1/year.

- effort:

  Fishing effort at which to assess the state.

- scope:

  `"whole"` or `"conditional"`.

- fixed:

  Owners held fixed for a conditional assessment. Use `"consumers"`,
  `"resource"`, or `"component:name"` for a component named `name`. This
  must be empty for whole-model assessments and nonempty for conditional
  assessments.

- criterion:

  `"state"` for the absolute state residual, or `"biomass"` for net
  biomass drift in the built-in checks.

- tolerances:

  Optional named numeric vector overriding individual tolerances. Names
  are `owner/measure`, as shown in the `measures` table.

- atol:

  Absolute tolerance for built-in checks, in grams per year.

- t:

  Time at which to evaluate the model, in years.

- dt:

  Step length made available to extension convergence callbacks.

- ...:

  Additional arguments made available to convergence callbacks.

## Value

A `MizerSteadyAssessment` list. `verdict` is `TRUE` when every required
owner is covered and passes, `FALSE` when any measured check fails, and
`NA` when no check fails but required support is missing. The list
records `scope`, `fixed`, `prescribed`, `required`, `covered`,
`missing`, `time` and `criterion`. `prescribed_masks` records the exact
prescribed coordinates for each owner. `unavailable` gives each missing
owner, capability and reason. `measures` has columns `owner`, `measure`,
`value`, `tolerance`, `units` and `passed`. A scope with no required
dynamic owners is steady by construction and has no measures.

## Details

The default `criterion = "state"` measures the integral of the absolute
rate of change of each consumer spectrum and the resource spectrum,
weighted by individual mass. Changes in different size classes therefore
cannot cancel. `criterion = "biomass"` instead takes the absolute value
after integrating and explicitly permits such cancellation. Both
built-in checks compare against `atol + tol * integral(w * abs(n))`.
They use the model's quadrature scheme and report values and tolerances
in grams per year.

Extensions supply named quantitative measures with their own units and
tolerances. Mizer compares `abs(value) <= tolerance` for every row and
preserves signed values in the report. `tol` and `atol` apply only to
the built-in checks. `tolerances` overrides measures from any provider.

## See also

[`isSteady()`](https://sizespectrum.org/mizer/reference/isSteady.md),
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md)
