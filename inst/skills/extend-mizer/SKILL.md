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

**Never let a rate jump as a function of abundance.** Depending on `t` or `w`
discontinuously is fine; depending on `n`, `n_pp` or `n_other` discontinuously
is not. Mizer's time steppers freeze the rates during each density update, so
they cannot see a threshold being crossed within a step. A rule like
`if (biomass < threshold) effort <- 0` gives a trajectory that keeps changing as
`dt` is refined, makes `steadyNewton()` stall, and makes `getStability()` return
a confident but meaningless answer — with no warning from mizer. Choosing
`method = "tr_bdf2"` does not help. Give the switch a finite width instead:

```r
# Bad: jumps. Good: ramps linearly between two thresholds.
frac <- (biomass - b_lim) / (b_trigger - b_lim)
effort <- effort * min(1, max(0, frac))
```

`max()`/`min()` kinks are continuous and much less troublesome — they cost some
accuracy, not correctness. See the
[Discontinuous rate functions](discontinuous_rates.html) article.

## Respecting the model's quadrature scheme

A model may be on either of two quadrature schemes, selected by the
`bin_average` entry of `second_order_w()` (see the `run-simulation` skill). Code
that ignores this looks correct and passes its tests, because the default is the
first-order scheme — and is then silently wrong by around 10% for anyone who has
switched the second-order scheme on. Three rules:

- **Build derived quantities from the rate functions, don't re-derive them.**
  `getEncounter()`, `getFeedingLevel()`, `getPredRate()`, `getEGrowth()` already
  carry the right quadrature for whichever scheme the model is on. Re-deriving a
  rate from the species parameters is how this goes wrong.
- **To go inside the encounter convolution, use `encounter_kernel()`, not
  `pred_kernel()`.** `pred_kernel()` returns the kernel point-sampled on the
  grid — right for plotting, and the form in which you *supply* a custom kernel,
  but not the bin-integrated coefficients the convolution consumes. Pair
  `encounter_kernel()` with the plain point prey weight
  `params@w_full * params@dw_full`, which is a normalisation the kernel is built
  to cancel, so it must not itself be bin-averaged.
- **If your setter precomputes an array from a size-dependent parameter,** gate
  it on `params@second_order_w[["bin_average"]]` the way `setExtMort()` and
  `setResource()` do, and do the integral once at setup so projection cost is
  unchanged.

For a summary-style integral $\int N_i(w) K_i(w) dw$, do not write the sum at
all: `sizeIntegral(params, weight = K, min_w = ..., max_w = ...)` does the
integral under whichever scheme the model is on and wraps the result — see
"Writing your own indicator" in the `analyse-and-plot` skill. If you need the
gating on its own, for a weight you are not integrating, that is
`bin_average_weight(K, params)`. Test any new integral with the flag **on** as
well as off; a test on the default path alone proves nothing.

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
- If a custom rate depends on the abundances through a threshold, read
  [Discontinuous rate functions](discontinuous_rates.html) before trusting any
  results.
