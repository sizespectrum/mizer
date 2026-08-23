---
name: extend-mizer
description: >-
  Extend or customise mizer's dynamics — add external food or mortality, replace
  a built-in rate calculation, or add an ecosystem component. Use for
  setExtEncounter(), setExtDiffusion() or setExtMort(); replacing
  mizerEncounter(), mizerPredRate(), mizerMort(), mizerEReproAndGrowth() or
  another rate with setRateFunction(); new dynamical pools with setComponent();
  extension subclasses; or second_order_w-aware custom quadrature. Pick the
  lightest mechanism that works. To change only an existing rate's parameters
  use the change-parameters skill; to load, save or share a model using an
  existing extension package use the use-extension-packages skill.
---

# Extending mizer

This is for changing how mizer works **without editing the package source**.
Pick the **lightest** mechanism that expresses your change, and reach for a
heavier one only when the lighter one cannot do the job. All setters return a new
`MizerParams` — reassign the result.

| Goal | Use |
|---|---|
| Add a fixed external food or mortality source (no new state variable) | `ext_encounter()` / `ext_mort()` |
| Change how one built-in rate is calculated | `setRateFunction()` |
| Add a new dynamical pool (detritus, carrion, oxygen, second resource…) | `setComponent()` |
| Store parameters for your custom code | `other_params()` (model-wide) or `component_params` (one component) |
| Extend plots/summaries for a custom model type | S3 methods on an S4 subclass |
| Replace arbitrary internal mizer code (last resort) | `customFunction()` |

If you only need to change one rate, prefer `setRateFunction()` over replacing
`mizerRates()` or patching internal functions. If you only need to change a
*number* the rates are computed from, you do not need this skill at all — see the
`change-parameters` skill.

## External encounter and mortality (lightest)

For a background process that affects fish but needs no state of its own. Both
are **species × size** arrays: `ext_encounter` (mass/year, added to
`getEncounter()`) and `ext_mort` (1/year, added to mortality).

```r
ext_mort(params) <- my_mort_array        # e.g. outside predators
ext_encounter(params) <- my_food_array   # extra unmodelled food
```

Build the array from the model's own grid rather than from literal dimensions,
so its shape and dimnames come out right, and **add to** what is already there
rather than overwriting it. An extra food source scaling allometrically with
body size:

```{r ext-encounter-example}
params_ext <- NS_params
extra_food <- outer(rep(0.1, nrow(species_params(params_ext))),
                    w(params_ext)^(3/4))
ext_encounter(params_ext) <- ext_encounter(params_ext) + extra_food
```

This is the right choice when the extra process is not itself depleted or
replenished, and the fish do not feed back on it. If it needs to respond to the
fish, you need a component instead.

## Replacing a rate function

Use when the model still follows mizer's standard flow but one step should be
computed differently — a time-dependent encounter, an alternative growth or
recruitment formulation, and so on. Mizer stores the **function name**, so the
function must be findable by name in the global environment or an installed
package; a function defined inside another function cannot be found.

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
`dt` is refined, makes `solver = "newton"` stall, and makes `getStability()` return
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

### Signatures and return shapes

Each replaceable rate has its own signature and its own expected return shape.
`setRateFunction()` calls your function with test inputs at registration and
checks the dimensions, so a mistake here surfaces immediately rather than
mid-projection.

| Rate | Signature | Return value |
|---|---|---|
| `Encounter` | `function(params, n, n_pp, n_other, t, ...)` | numeric matrix, species × size |
| `FeedingLevel` | `function(params, n, n_pp, n_other, t, encounter, ...)` | numeric matrix, species × size |
| `EReproAndGrowth` | `function(params, n, n_pp, n_other, t, encounter, feeding_level, ...)` | numeric matrix, species × size |
| `ERepro` | `function(params, n, n_pp, n_other, t, e, ...)` | numeric matrix, species × size |
| `EGrowth` | `function(params, n, n_pp, n_other, t, e_repro, e, ...)` | numeric matrix, species × size |
| `PredRate` | `function(params, n, n_pp, n_other, t, feeding_level, ...)` | numeric matrix, species × **full** size grid |
| `PredMort` | `function(params, n, n_pp, n_other, t, pred_rate, ...)` | numeric matrix, species × size |
| `FMort` | `function(params, n, n_pp, n_other, t, effort, e_growth, pred_mort, ...)` | numeric matrix, species × size |
| `Mort` | `function(params, n, n_pp, n_other, t, f_mort, pred_mort, ...)` | numeric matrix, species × size |
| `RDI` | `function(params, n, n_pp, n_other, t, e_growth, mort, e_repro, diffusion, ...)` | numeric vector, one value per species |
| `RDD` | `function(rdi, species_params, params, t, ...)` | numeric vector, one value per species |
| `ResourceMort` | `function(params, n, n_pp, n_other, t, pred_rate, ...)` | numeric vector, one value per full size bin |
| `Diffusion` | `function(params, n, n_pp, n_other, t, feeding_level, ...)` | numeric matrix, species × size |
| `Rates` | `function(params, n, n_pp, n_other, t, effort, rates_fns, ...)` | named list with all standard rate components |

Three rules that follow from the table:

- **Return plain numeric objects.** Never return an `ArraySpeciesBySize` or
  `ArrayTimeBySpecies` from a function registered with `setRateFunction()`; the
  `get…()` wrappers add those classes afterwards where they belong.
- **Most size-resolved rates match `initialN(params)`** in both dimensions and
  dimnames. The exceptions are `PredRate`, which is on the full prey grid
  (`w_full`), and `ResourceMort`, which is one value per `w_full` bin rather than
  one per species.
- **Build outputs from existing mizer arrays** so dimensions and dimnames are
  inherited. Most bugs in extension code are the right numbers in the wrong
  shape.

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
all: `sizeIntegral(params, weighting = K, min_w = ..., max_w = ...)` does the
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

The functions named in that call take the component's name as a `component`
argument, so one implementation can serve several components, and reach their
own state as `n_other[[component]]` and their own parameters as
`params@other_params[[component]]`. A dynamics function additionally receives
the current `rates` and the step length `dt`, and returns the component's **new
state** — not a rate of change, so integrate the step yourself. This pair makes
a detritus pool that fish eat and that relaxes back towards a capacity:

```{r detritus-functions}
detritusEncounter <- function(params, n, n_pp, n_other, component, ...) {
    params2 <- params
    # Drop this component before delegating, or mizerEncounter() calls back
    # into this function and recurses.
    params2@other_encounter[[component]] <- NULL
    mizerEncounter(params2, n = n, n_pp = n_other[[component]],
                   n_other = n_other, ...)
}

detritusDynamics <- function(params, n_other, rates, dt, component, ...) {
    detritus <- n_other[[component]]
    p <- params@other_params[[component]]
    interaction <- params@species_params$interaction_resource
    mort <- as.vector(interaction %*% rates$pred_rate)
    target <- p$rate * p$capacity / (p$rate + mort)
    # Exact over the step, so the result does not depend on dt
    target - (target - detritus) * exp(-(p$rate + mort) * dt)
}
```

Passing the component's own abundance to `mizerEncounter()` as `n_pp`, as
`detritusEncounter()` does, is the trick for making a component act as an extra
prey spectrum: it is then eaten through the ordinary predation kernel and shows
up in `getDiet()` without further work.

<!-- article-only -->

## Worked example: external encounter and mortality

The `extra_food` built above adds straight to the total encounter rate:

```{r ext-encounter-check}
enc_base <- getEncounter(NS_params)
enc_ext <- getEncounter(params_ext)
range(enc_ext - enc_base, na.rm = TRUE)
```

External mortality works the same way — note the negative exponent, since
mortality falls with size where the extra food rose with it:

```{r ext-mort-example}
params_mort <- NS_params
extra_mort <- outer(rep(0.05, nrow(species_params(params_mort))),
                    w(params_mort)^(-1/4))
ext_mort(params_mort) <- ext_mort(params_mort) + extra_mort
```

## Worked example: a seasonal encounter rate

Wrapping a built-in rate, in full. This one takes the amplitude and period of a
seasonal cycle from `other_params()`:

```{r seasonal-encounter-setup}
params <- NS_params
other_params(params) <- list(season_amplitude = 0.2, season_period = 1)
```

```{r seasonal-encounter-fun}
seasonalEncounter <- function(params, n, n_pp, n_other, t, ...) {
    p <- other_params(params)
    multiplier <- 1 + p$season_amplitude * sin(2 * pi * t / p$season_period)
    multiplier * mizerEncounter(params, n = n, n_pp = n_pp, n_other = n_other,
                                t = t, ...)
}
```

Registered by name, the rate now moves with `t`:

```{r seasonal-encounter-register}
params2 <- setRateFunction(params, "Encounter", "seasonalEncounter")
enc0 <- getEncounter(params2, t = 0)
enc_quarter <- getEncounter(params2, t = 0.25)
range(enc_quarter / enc0, na.rm = TRUE)
```

At `t = 0.25` the multiplier is at its maximum, `1 + season_amplitude`, exactly
as intended. The seasonality then carries through a projection:

```{r seasonal-biomass}
sim <- project(params2, t_max = 2, t_save = 0.1)
plotBiomass(sim)
```

A good first test for such a function checks that at `t = 0` it agrees with
`mizerEncounter()`, that at `t = 0.25` it is scaled by `1 + season_amplitude`,
and that its result has the same dimensions and dimnames as `initialN(params)`.

## Worked example: a detritus-like component

Putting `detritusEncounter()` and `detritusDynamics()` from above to work. The
component is stored on the full resource size grid, and starts at half the
capacity it relaxes towards:

```{r detritus-component}
detritus_params <- list(capacity = initialNResource(params),
                        rate = params@rr_pp)

params3 <- setComponent(
    params,
    component = "Detritus",
    initial_value = initialNResource(params) / 2,
    dynamics_fun = "detritusDynamics",
    encounter_fun = "detritusEncounter",
    component_params = detritus_params,
    colour = "orange"
)
```

Its initial state is now in `initialNOther(params3)$Detritus` and its settings in
`getComponent(params3, "Detritus")`. Because it is eaten through the ordinary
predation kernel, it appears in the diet with no further work:

```{r detritus-diet}
plotDiet(params3, species = "Cod")
```

To let the component kill fish as well as feed them, give `setComponent()` a
`mort_fun` too.

<!-- /article-only -->

## Extending plots and summaries: S3 methods on an S4 subclass

This is the most flexible route that still works with mizer's public generics.
Although `MizerParams` and `MizerSim` are S4 classes, mizer registers most of its
user-facing methods as **S3** methods. So an extension defines a formal S4
subclass of these objects and provides S3 methods for it —
`getBiomass.MyMizerSim()`, `plotBiomass.MyMizerSim()`,
`summary.MyMizerParams()` — which is what makes a summary or plot account for
components the extension added. Every such method must call `NextMethod()`, so
that several extensions loaded at once compose instead of overwriting each other.

This only really works inside a package, because the marker class is created by
mizer when the package loads. Writing that package — the class, the `.onLoad`
registration and the methods — is the subject of the
`create-extension-package` skill, and the order the methods then run in is the
subject of the `use-extension-packages` skill.

## Storing parameters

- **`other_params(params)$foo <- value`** — model-wide parameters your custom
  functions can read.
- **`component_params`** (a `setComponent()` argument) — parameters scoped to one
  component.

Keeping the two separate is what keeps `params@other_params` readable when
several custom functions are in play.

## `customFunction()` (last resort)

`customFunction()` replaces an internal mizer function inside the package
namespace. Reach for it only after confirming that `setRateFunction()`,
`setComponent()`, `resource_dynamics()<-` and `setReproduction()` cannot express
the change: a replacement that is not fully compatible breaks the package, and
nothing checks that it is.

## Testing and packaging

- Test **both** halves: that `setRateFunction()` or `setComponent()` succeeds,
  and that the downstream function you actually care about — `getEncounter()`,
  `getDiet()`, `project()` — uses your extension as intended.
- After registering a custom rate, run `project()` on a short horizon and check
  the affected `get…()` output looks right; write a `testthat` test that compares
  against a hand-computed value on a small model.
- If a custom rate depends on the abundances through a threshold, read
  [Discontinuous rate functions](discontinuous_rates.html) before trusting any
  results.
- Once an extension is useful in more than one project, or you want to give it
  to someone else, make it a package: a stable namespace mizer can resolve
  function names in, somewhere for tests and documentation to live, and a version
  that gets recorded in `params@extensions`. Everything that only matters once
  you share — `.onLoad` registration, marker classes, dispatch via
  `NextMethod()`, bundled data objects, reporting through `info_level`, and
  upgrading objects saved by an earlier version — is in the
  `create-extension-package` skill. For the other side of it, loading and saving
  a model that needs an extension package, see the
  `use-extension-packages` skill.
- For a larger design, it is worth discussing the interface on the
  [mizer issue tracker](https://github.com/sizespectrum/mizer/issues) before
  committing to it.
