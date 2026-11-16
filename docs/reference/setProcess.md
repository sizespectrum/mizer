# Declare capabilities for a model process

Register a convergence callback or a prescribed-input declaration for
the current consumers or resource process. The declaration is stored in
the model together with the identity of its current dynamics provider.
Replacing that provider or changing the projection dispatch chain
invalidates the declaration until support is explicitly adopted again.
Changing parameters or consumer rate functions leaves the identity of
the fish equations intact. Calling `setProcess(params, process)` with no
capability arguments explicitly adopts the current provider for all
retained capability families. A partial patch adopts only the families
explicitly supplied: convergence/role, stationary/equation rates/tuning,
reproduction level, transformations, or mortality diagnostics,
prescribed inputs, temporal dependence or resource equations. Omitted
families keep their previous provider bindings, including staleness. A
replacement definition adopts all its supplied capabilities.

## Usage

``` r
setProcess(
  params,
  process,
  definition,
  convergence_fun,
  role,
  stationary,
  equation_rates_fun,
  tuning,
  level_fun,
  transformations,
  mortality_fun,
  prescribed_mask_fun,
  autonomous,
  rate_fun
)

getProcess(params, process)
```

## Arguments

- params:

  A MizerParams object.

- process:

  `"consumers"`, `"resource"` or `"reproduction"`. The latter has no
  state role or convergence callback.

- definition:

  A named list of assessment and tuning capabilities.

- convergence_fun:

  Name of an assessment callback, optionally qualified as
  `package::exported_function`, or `NULL`.

- role:

  Either `"dynamic"` or `"prescribed"` for a state process.

- stationary:

  Explicit TRUE/FALSE stationary-tuning declaration.

- equation_rates_fun:

  Callback exposing effective equation quantities.

- tuning:

  Named strategies and their default, or NULL to remove support.

- level_fun:

  Law-specific reproduction-level callback, or NULL.

- transformations:

  Named operation declarations for structural and scale changes. Each
  declares action = 'none' or a named transform_fun with action =
  'update', explicit paths and optional settings. Supply this argument
  or a new definition to adopt support for a changed provider; an
  assessment-only update retains the previous transformation binding.

- mortality_fun:

  Named consumer mortality-breakdown callback, or NULL to remove
  diagnostic support. It is bound independently of other families.

- prescribed_mask_fun:

  Named structural mask callback. It receives params and optional
  owner/process and returns logical values with the exact state
  dimensions and names; TRUE marks prescribed inputs. It is not called
  with time or state values. Use role = 'dynamic' for a partial mask.

- autonomous:

  TRUE for equations without explicit time dependence, FALSE for forced
  equations, or NULL for undeclared support. This is independent of
  stationary tuning. Declare reproduction autonomy on consumers, whose
  temporal binding includes the complete effective rate pipeline.

- rate_fun:

  Named exact resource equation callback returning dn_pp/dt, with the
  state dimensions and names. Its context includes params, n, n_pp,
  n_other, rates and t. Only the resource process accepts this
  capability.

## Value

`setProcess()` returns the updated model. `getProcess()` returns the
callback, role, reusable `definition` and recorded provider `binding`.
With no `process` argument, `getProcess()` returns all stored processes.

## Details

A convergence callback returns a data frame with the columns `measure`,
`value`, `tolerance` and `units`. It supplies measurements, not
verdicts;
[`getSteadyAssessment()`](https://sizespectrum.org/mizer/reference/getSteadyAssessment.md)
compares each absolute value with its tolerance. Callback arguments
include the model, state, complete rate list, derivative, time and
process name; accept `...` for arguments not needed by the callback.

A new `definition` replaces the previous callback and role before
explicit arguments are applied. Without `definition`, omitted arguments
are retained. The default role is `"dynamic"`. A `"prescribed"` role
declares the whole state block as an input rather than a dynamic
variable. A dynamic owner may declare a partial `prescribed_mask_fun`.
Supplying `convergence_fun = NULL` removes its callback. State and
equations are not changed by this registration.

## See also

[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`getSteadyAssessment()`](https://sizespectrum.org/mizer/reference/getSteadyAssessment.md)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`mizerReproductionDefinition()`](https://sizespectrum.org/mizer/reference/mizerReproductionDefinition.md),
[`other_mort()`](https://sizespectrum.org/mizer/reference/other_mort.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
