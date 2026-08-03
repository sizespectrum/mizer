---
name: extend-mizer
description: >-
  Extend or customise mizer's dynamics — add external food/mortality, replace a
  built-in rate calculation, or add a new ecosystem component. Use whenever the
  user wants a custom encounter/growth/mortality/reproduction formulation
  (setRateFunction), a background food or predation source (setExtEncounter,
  setExtMort), a new dynamical pool like detritus or carrion (setComponent), or
  asks how to make mizer do something its standard setters do not cover.
---

# Extending mizer

Pick the **lightest** mechanism that expresses your change. Reach for a heavier
one only when the lighter one cannot do the job. All setters return a new
`MizerParams` — reassign the result.

| Goal | Use |
|---|---|
| Add a fixed external food or mortality source (no new state variable) | `setExtEncounter()` / `setExtMort()` |
| Change how one built-in rate is calculated | `setRateFunction()` |
| Add a new dynamical pool (detritus, carrion, oxygen, second resource…) | `setComponent()` |
| Store parameters for your custom code | `other_params()` (model-wide) or `component_params` (one component) |
| Extend plots/summaries for a custom model type | S3 subclass + new methods |
| Replace arbitrary internal mizer code (last resort) | `customFunction()` |

## External encounter and mortality (lightest)

For a background process that affects fish but needs no state of its own. Both
are **species × size** arrays: `ext_encounter` (mass/year, added to
`getEncounter()`) and `ext_mort` (1/year, added to mortality).

```r
params <- setExtMort(params, ext_mort = my_mort_array)            # e.g. outside predators
params <- setExtEncounter(params, ext_encounter = my_food_array)  # extra unmodelled food
```

## Replacing a rate function

Use when the model still follows mizer's standard flow but one step should be
computed differently — a time-dependent encounter, an alternative growth or
recruitment formulation, and so on. Mizer stores the **function name**, so the
function must be findable by name in the global environment or an installed
package.

```r
params <- setRateFunction(params, "Mort", "myMort")
getRateFunction(params)      # list the current rate functions
```

Replaceable rates include `Rates`, `Encounter`, `FeedingLevel`, `PredRate`,
`PredMort`, `Mort`, `ResourceMort`, `EReproAndGrowth`, `ERepro`, `EGrowth`,
`Diffusion`, `FMort`, `RDI` (density-independent recruitment), and `RDD`
(density-dependent recruitment). Resource dynamics is set separately via
`resource_dynamics()`.

**Write your function by starting from the built-in** (`mizerMort()`,
`mizerEncounter()`, …): copy it, change what you need, and keep the same
signature and return dimensions/dimnames. Custom rate functions receive
`params`, the current state (`n`, `n_pp`, `n_other`), `t`, and previously
computed rates via `...` — accept `...` and pull what you need from it.

```r
# Add a size-independent extra mortality, its size stored in other_params
myMort <- function(params, n, n_pp, n_other, t, f_mort, pred_mort, ...) {
    base <- mizerMort(params, n, n_pp, n_other, t, f_mort, pred_mort, ...)
    base + other_params(params)$extra_mort  # a scalar or species × size array
}
other_params(params)$extra_mort <- 0.1     # store the parameter
params <- setRateFunction(params, "Mort", "myMort")
```

**Time-dependent rates are the key reason to reach for `setRateFunction()`.**
Species parameters and rate arrays are fixed for the whole simulation, but a rate
*function* receives the current time `t` and can therefore change as the run
proceeds — seasonal forcing, a warming trend, a management measure that switches
on in a given year. Wrap the built-in and scale its result by `t`:

```r
seasonalMort <- function(params, t, ...) {
    mizerMort(params, t = t, ...) * (1 + 0.3 * sin(2 * pi * t))   # t in years
}
params <- setRateFunction(params, "Mort", "seasonalMort")
```

## Adding a component

Use `setComponent()` for a new dynamical quantity — any R object: a scalar, a
vector on the size grid, or a list. A component may contribute in up to three
ways:

- `dynamics_fun` — updates the component's own state during projection,
- `encounter_fun` — adds to `getEncounter()`,
- `mort_fun` — adds to `getMort()`.

```r
params <- setComponent(
    params, "detritus",
    initial_value    = 1e5,
    dynamics_fun     = "detritus_dynamics",
    encounter_fun    = "detritus_encounter",  # optional
    component_params = list(rho = 0.1)
)
```

Access components with `getComponent()`, remove them with `removeComponent()`.
Component state is available to all custom functions as `n_other[["detritus"]]`,
and its parameters via `component_params`.

## Storing parameters

- **`other_params(params)$foo <- value`** — model-wide parameters your custom
  functions can read.
- **`component_params`** (a `setComponent()` argument) — parameters scoped to one
  component.

## Testing and packaging

- After registering a custom rate, run `project()` on a short horizon and check
  the affected `get…()` output looks right; write a `testthat` test that compares
  against a hand-computed value on a small model.
- To share an extension as an R package (with `.onLoad` registration and method
  dispatch via `NextMethod()`), follow the
  [Creating extension packages](creating-extension-packages.html) article.
- For the full treatment of everything above, see the
  [Extending mizer](extending-mizer.html) article.
