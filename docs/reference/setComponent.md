# Add a dynamical ecosystem component

By default, mizer models any number of size-resolved consumer species
and a single size-resolved resource spectrum. Your model may require
additional components, like for example detritus or carrion or multiple
resources or .... This function allows you to set up such components.

## Usage

``` r
setComponent(
  params,
  component,
  initial_value,
  dynamics_fun,
  encounter_fun,
  mort_fun,
  component_params,
  colour = "grey",
  linetype = "solid",
  definition,
  rate_fun,
  integration_order,
  convergence_fun,
  role,
  stationary,
  equation_rates_fun,
  tuning,
  transformations,
  biomass_fun,
  biomass_filters,
  prescribed_mask_fun,
  autonomous
)

removeComponent(params, component)

getComponent(params, component)
```

## Arguments

- params:

  A MizerParams object

- component:

  Name of the component of interest. If missing, a list of all
  components will be returned.

- initial_value:

  Initial value of the component. Required for a new component; retained
  when omitted for an existing component. Must not be `NULL`; use
  `removeComponent()` to remove a component.

- dynamics_fun:

  Name of function to calculate value at the next time step. Required
  unless `rate_fun` is supplied.

- encounter_fun:

  Name of function to calculate contribution to encounter rate.
  Optional.

- mort_fun:

  Name of function to calculate contribution to the mortality rate.
  Optional.

- component_params:

  Object holding the parameters needed by the component functions. This
  could for example be a named list of parameters. Optional.

- colour:

  Line colour to use for the component in plots. Defaults to `"grey"`.

- linetype:

  Line type to use for the component in plots. Defaults to `"solid"`.

- definition:

  A named list specifying the component's callbacks and plotting styles.
  Supplying a list replaces its previous definition before explicitly
  supplied arguments are applied. State and component parameters are
  supplied separately. Updating only state or component parameters does
  not re-adopt capabilities after a global provider change; supply the
  relevant capability explicitly.

- rate_fun:

  Name of function returning the component's derivative. Required unless
  `dynamics_fun` is supplied.

- integration_order:

  Optional positive integer declaring the order of a custom
  `dynamics_fun`. `NULL` removes the declaration.

- convergence_fun:

  Name of a convergence-measurement callback, or `NULL` to remove it. A
  dynamic component without this support makes a whole-model assessment
  incomplete.

- role:

  Either `"dynamic"` (the default) or `"prescribed"`.

- stationary:

  Explicit TRUE/FALSE stationary-tuning declaration.

- equation_rates_fun:

  Callback exposing effective equation quantities.

- tuning:

  Named strategies and their default, or NULL to remove support.

- transformations:

  Named transformation declarations. Update hooks may modify this
  component's initial state and parameters; additional paths must be
  explicit. Use action = 'none' for supported no-action operations.

- biomass_fun:

  Named biomass measurement callback, or NULL to remove measurement
  support. This is independent of convergence and dynamics.

- biomass_filters:

  Character vector of explicitly supported filters: min_w, max_w, min_l,
  max_l and use_cutoff. Defaults to none.

- prescribed_mask_fun:

  Named structural input-mask callback, or NULL. For a dynamic numeric
  component, return logical values with exactly its state dimensions and
  names; TRUE marks prescribed inputs. The callback receives params and
  optional component/owner, without time or state.

- autonomous:

  TRUE for equations and contributions without explicit time dependence,
  FALSE for forced providers, or NULL for undeclared support. This
  declaration is independent of stationary tuning and convergence.

## Value

The updated MizerParams object

For `getComponent`: A list with the entries `initial_value`,
`dynamics_fun`, `rate_fun`, `encounter_fun`, `mort_fun`,
`component_params`, `convergence_fun`, `role`, and `definition` for the
requested component. If the requested component does not exist, `NULL`
is returned. If no `component` argument is given, then a list of lists
for all components is returned.

## Details

The component can be a number, a vector, an array, a list, or any other
data structure you like.

A new component needs an `initial_value`, supplied separately from its
definition. Supply at least one of `dynamics_fun` and `rate_fun`,
directly or through `definition`. For an existing component, omitted
arguments retain their existing values. Explicit `NULL` removes an
optional callback, provided at least one of the two dynamics
capabilities remains. Remove a component with `removeComponent()`.

A `definition` is a reusable named list containing `dynamics_fun`,
`rate_fun`, `encounter_fun`, `mort_fun`, `colour`, `linetype` and/or
`integration_order`, `convergence_fun`, `role`, `stationary`,
`equation_rates_fun` and/or `tuning`. Explicit arguments override its
entries. Supplying a definition replaces the component's previous
definition: omitted optional callbacks are removed and omitted styles
default to `"grey"` and `"solid"`. Without `definition`, explicitly
supplied arguments patch the existing definition instead. Initial state
and `component_params` belong to the component instance and are retained
when omitted, including when its definition is replaced.

The assembled definition is stored in the model. Subsequent changes to
the list used as `definition` do not change that model; call
`setComponent()` again to apply an update. Callback function names are
stored, rather than their implementations, so the functions still need
to be available when the model is used.

## Equations and custom steps

A `rate_fun` returns the component's derivative, with the numeric shape
and names of its state preserved. It receives `params`, `n`, `n_pp`,
`n_other`, `rates`, `t`, `component` and `...`. Mizer applies Euler
under first-order projection and Heun under the second-order methods.
This is explicit integration, including when the fish use TR-BDF2.

A `dynamics_fun` returns the new state and owns the integration. It is
the appropriate capability for stiff or non-numeric states. When both
callbacks are supplied, the custom step performs projection and the
derivative remains available for diagnostics. Callback names may be
unqualified or use `package::exported_function`.

Custom steps may accept `predicted_state` and `endpoint_rates`. During a
second-order corrector these contain the predicted endpoint state (a
list with `n`, `n_pp`, `n_other` and `t`) and its rate list. They are
`NULL` when no prediction is available. The ordinary state arguments
still describe the start of the step. Prediction data do not by
themselves establish the accuracy or stability of the custom
integration.

`integration_order` declares the order of a custom step; omission means
unknown, not an assumed second order. Replacing `dynamics_fun` clears an
inherited order unless an order is explicitly supplied with the
replacement. Mizer chooses the integrator for equations-only components,
so they must not declare a custom integration order.

## Steady-state assessment

A `convergence_fun` returns a data frame with `measure`, `value`,
`tolerance` and `units` columns. It receives the model, state, complete
rate list, time, derivative when available and component name. Accept
`...` for unused arguments. Mizer compares absolute measured values to
their tolerances; the callback does not return its own pass/fail
verdict.

Components are dynamic by default, including legacy components using
`constant_other`. Declare `role = "prescribed"` explicitly when the
entire component is a prescribed input. Temporarily holding a dynamic
component fixed during analysis requires conditional assessment scope.
Replacing a dynamics callback removes its inherited convergence callback
and resets an inherited prescribed role to dynamic unless those
declarations are explicitly supplied again with the replacement.

## See also

"Extending mizer": [guide to extending
mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.html)

Other extension tools:
[`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md),
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md),
[`initialNOther<-()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md),
[`mizerReproductionDefinition()`](https://sizespectrum.org/mizer/reference/mizerReproductionDefinition.md),
[`other_mort()`](https://sizespectrum.org/mizer/reference/other_mort.md),
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md),
[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md),
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
