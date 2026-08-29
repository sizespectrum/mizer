---
name: upgrade-extension-package
description: >-
  Bring an existing mizer extension package up to date with a newer mizer. Use
  when a package written against an earlier mizer needs migrating to the S3
  MizerParams and the simplified extension mechanism; when its methods have
  stopped dispatching, its S4 marker classes or registerExtension() calls have
  to go, or its tests fail on mizer's new reports; or when a workaround it
  carries is no longer needed. Only for a genuine extension — one that declares
  marker classes, registers methods on mizer's generics, installs a rate
  function or registers a component — and only in addition to the
  upgrade-mizer-code skill, which covers everything a package does as ordinary
  user code, including the version floor, docs and test suite of any package
  that merely calls mizer. For writing a new extension package use the
  create-extension-package skill.
---

# Upgrading your extension package

This article is for the maintainer of an existing mizer extension package —
`mizerShelf`, `mizerMR`, `mizerStarvation`, `therMizer` or your own — bringing it
up to date with a newer mizer. It has two halves:

- **Part 1** is the one-off migration to mizer's S3 architecture and simplified
  extension mechanism (mizer 3.4). Every package written against mizer 3.3 or
  earlier needs it.
- **Part 2** collects the other mizer changes that land specifically on
  extension code — methods that no longer dispatch, workarounds that can be
  dropped, tests that now fail. These are release-by-release and each one is
  independent of the S3 migration.

It covers only what is specific to *being an extension*.

## Before you start: is this an extension?

Not every package that depends on mizer is a mizer extension, and only an
extension needs this article. Look for any of these:

- `setClass()` declaring a marker class, or `.onLoad()` calling
  `registerExtension()` / `registerExtensions()`
- methods on mizer's own generics — `project.myClass`,
  `setParams.myClass`, `steadySingleSpecies.myClass`
- a rate function installed with `setRateFunction()`, or a component
  registered with `setComponent()`

A package with none of them — one that only *calls* mizer, however extensively,
like `mizerExperimental` — is ordinary user code in a package wrapper. The
`upgrade-mizer-code` skill covers it in full, including its "If your code is a
package" section. Stop here and use that.

**Everything in this article is in addition to that skill, not instead of it.**
Work through the `upgrade-mizer-code` skill first — its code-pattern index for
the asymptomatic changes, its symptom index for anything that broke, and its
package section for the version floor, the changelog, the docs and the test
suite. Then come back here for what only an extension has to do. For what the
finished package should look like, rather than how to get there, read the
`create-extension-package` skill: this article does not repeat its rules, and
its checklist for package authors is the acceptance test for a migrated
package.

<!-- agent-only -->

## Symptom index

| Symptom | Cause | Section |
|---|---|---|
| `is()` / `setClass()` / `expect_s4_class()` in the package, or `.onLoad()` calling `registerExtension()` / `getRegisteredExtensions()` | written against the S4 mizer | Part 1 |
| After `rm(list = ls())`, or in the second and later examples of `R CMD check`, mizer reports the extension's class is not a defined class | the dynamic S4 marker classes lived in `.GlobalEnv` | Part 1 — the mechanism is gone, migrate rather than repair |
| After `devtools::load_all()`, `coerceToExtensionClass()` reports no method or default for coercing `MizerParams` | same obsolete marker-class machinery | Part 1 |
| The package's results are unchanged when the model is switched to `second_order_w = TRUE`, or a diagnostic is out by a constant factor | the extension's own integrals ignore the `bin_average` flag | Working under both numerical schemes |
| A component's state advances twice as far per step, or its `dynamics_fun` double-counts, under `method = "second_order"` (`"tr_bdf2"`) or `"predictor_corrector"` | the second-order steppers call `dynamics_fun` twice per step, from the same start-of-step state | Working under both numerical schemes |
| A method the package defines is never called, and its name starts with `get` | the `get`-prefixed accessors are plain aliases now; dispatch happens on the bare name | Methods on the renamed rate-array accessors |
| `steady.MyClass` runs but `tuneSteadyState()` ignores it | the new names are separate generics | Methods on the steady-state finders |
| The package removes a species parameter column by writing into `params@species_params` | supported route exists since 3.3.1 | Withdrawing a species parameter column |
| The package declares the current `gamma` as a given species parameter to stop it moving | the default is measured on mizer's own reference state now | The `gamma`/`f0` workaround can go |
| Tests fail on message wording, or on a call that used to be silent | new reports and rephrased messages | Tests that match messages or assert silence |
| Links to `articles/creating-extension-packages.html`, `articles/extending-mizer.html`, `vignette("cheatsheet-…")` | the articles were renamed to guides | Links to renamed articles and skills |

<!-- /agent-only -->

## Part 1 — the S3 extension mechanism

### What changed in mizer

Starting with mizer 3.3.0.9000 / 3.4:

1. **`MizerParams` and `MizerSim` are S3 classes.** They are stored as lists with
   class vectors `c("MizerParams")` and `c("MizerSim")`.
2. **S3 extension classes:** an extension package prepends its class name to the
   S3 class vector — `c("mizerShelf", "MizerParams")` for params, and
   `c("mizerShelfSim", "MizerSim")` for simulations.
3. **Session registration is obsolete:** `registerExtension()` and the dynamic S4
   marker classes created at package load time have been removed.
4. **No `.onLoad` hook or active bindings needed:** bundled datasets saved as S3
   objects retain their class vectors and extension metadata under standard R
   lazy-loading.
5. **Extension tracking via `recordExtension()`:** the extension requirement and
   package version are recorded in `params$extensions` (or `params@extensions`).
6. **Class coercion via `coerceToExtensionClass()`:** it reads `params$extensions`
   and sets the S3 class vector accordingly. `project()` automatically coerces
   its `MizerSim` output to the corresponding `…Sim` extension class.

### Step by step

#### 1. Clean up `R/<pkgname>-package.R`

Remove `.onLoad()` entirely (and delete `R/zzz.R` if it held only `.onLoad()`), or the `registerExtension()` and
`makeActiveBinding()` calls within it if it does anything else.

- **Remove `registerExtension()`**: session registration is no longer used.
- **Remove `makeActiveBinding()`**: datasets in `data/` now carry their S3 class
  directly.
- **Remove `@importFrom methods is`**: replace `is(x, "MizerParams")` in code
  with `inherits(x, "MizerParams")`.
- **Remove deleted files from `Collate:`**: if `DESCRIPTION` has a manual `Collate:` field, remove `zzz.R` (or whichever file held `.onLoad()`) from it.

```r
# BEFORE:
.onLoad <- function(libname, pkgname) {
    mizer::registerExtension(pkgname, requirement = "sizespectrum/myExtension")
    if (exists("my_params", envir = asNamespace(pkgname), inherits = FALSE)) {
        ns <- asNamespace(pkgname)
        raw <- get("my_params", envir = ns)
        makeActiveBinding("my_params",
                          fun = function() mizer::coerceToExtensionClass(raw),
                          env = ns)
    }
}

# AFTER:
# Delete .onLoad() completely if there are no other load actions.
```

#### 2. Delete every marker class declaration

An extension class is now an ordinary entry in the object's S3 class vector, so
nothing has to be declared. Delete

```r
setClass("myExtension", contains = "MizerParams")
setClass("myExtensionSim", contains = "MizerSim")
```

and, if the package carried the workaround for the marker classes being wiped by
`rm(list = ls())` or by the `cleanEx()` that `R CMD check` runs before every
example, delete that too:

```r
setClass("myExtension", contains = "MizerParams", where = globalenv())
```

Both are actively harmful now: a sealed S4 class of the same name cannot be
re-parented, which is what prevented two independently developed extensions from
being chained in either load order.

Then check `NAMESPACE` and the roxygen comments for `@exportClass`,
`@importClassesFrom mizer` and `setValidity()` on those names, and remove those
as well. `devtools::document()` will not drop a stale `exportClasses()` line
while the `setClass()` is still there.

#### 3. Update the class documentation

Rewrite the roxygen block that documented the extension class:

```r
#' myExtension extension classes
#'
#' S3 extension classes for [MizerParams] and [MizerSim] that enable S3 dispatch
#' for extension-specific methods.
#'
#' The class names are ordinary entries in the object's S3 class vector. All
#' extension-specific data lives in `other_params(params)` or in component
#' parameters (see [setComponent()]).
#'
#' Objects of class `myExtension` are created by [newMyExtensionParams()].
#' Objects of class `myExtensionSim` are returned automatically by [project()]
#' when called on a `myExtension` params object.
#'
#' No class declaration is needed. [newMyExtensionParams()] records the
#' extension on the object with [mizer::recordExtension()] and then calls
#' [mizer::coerceToExtensionClass()].
#'
#' @name myExtension-class
#' @keywords internal
NULL
```

#### 4. Update constructors and setup functions

Every constructor and setup function should

1. call `recordExtension()` with the package name, `version` and explicit
   `requirement = "<owner>/<repo>"` (now required since session registration
   no longer supplies it), and
2. call `coerceToExtensionClass(params)`, for a dispatching extension.

```r
params <- mizer::recordExtension(
    params, "myExtension",
    version = as.character(utils::packageVersion("myExtension")),
    requirement = "sizespectrum/myExtension"
)
params <- mizer::coerceToExtensionClass(params)
```

For a **metadata-only** extension such as `mizerStarvation`, call
`recordExtension()` and omit `coerceToExtensionClass()`.

#### 5. Migrate bundled data objects in `data/`

`MizerParams` and `MizerSim` objects in `data/` that were stored as S4 objects
have to be upgraded to S3 and re-saved. Run this once, per object:

```r
library(mizer)
load("data/my_params.rda")

# Convert S4 to S3 list if needed
p <- mizer:::upgrade_s4_to_s3(my_params)
class(p) <- c("myExtension", "MizerParams")
p <- recordExtension(
    p, "myExtension",
    version = as.character(utils::packageVersion("myExtension")),
    requirement = "sizespectrum/myExtension"
)
p <- validParams(p)

my_params <- p
save(my_params, file = "data/my_params.rda", compress = "bzip2")
```

and update the documentation in `R/data.R`:

```r
#' @format A [myExtension-class] object.
#'
#' The object was created with [newMyExtensionParams()] and is stored in
#' the package's `data/` directory. R's standard lazy-loading preserves its S3
#' class vector and extension metadata, so no load hook or active binding is
#' needed.
```

#### 6. Replace S4 assertions and type checks

- Replace `is(x, "MizerParams")` and `is(x, "MizerSim")` with
  `inherits(x, "MizerParams")` and `inherits(x, "MizerSim")`.
- In `tests/testthat/`, replace `expect_s4_class()` with `expect_s3_class()`
  throughout — for `MizerParams` and `MizerSim` as well as for your own
  `myExtension` and `myExtensionSim`.
- Replace any assertions on `mizer::getRegisteredExtensions()` with assertions on
  `mizer::getMetadata(params)$extensions` (e.g., checking
  `getMetadata(params)$extensions$myExtension$version` and `$requirement`).

#### 7. Provide an `upgrade` method if layout has evolved

If your extension's internal layout has changed across versions, register an
idempotent S3 method for `utils::upgrade()` so that older saved objects are
migrated automatically when validated. See "Upgrading objects across versions of
your extension" in the `create-extension-package` skill:

```r
#' @exportS3Method utils::upgrade
upgrade.myExtension <- function(object, ...) {
    # Detect old layout and migrate idempotently
    object
}
```

#### 8. Update `DESCRIPTION`

- Remove `methods` from `Imports:` if it was there only for the S4 checks.
- Require `mizer (>= 3.3.0.9000)`, or `mizer (>= 3.4.0)` once that is released.
- If `DESCRIPTION` has an explicit `Collate:` field, update it to remove deleted
  files (such as `zzz.R`) and include any new ones (such as `upgrade.R`).

#### 9. Update vignettes, narrative docs and `NEWS.md`

- Wherever a vignette describes S4 marker classes, `.onLoad()` or
  `registerExtension()`, replace it with S3 extension classes,
  `recordExtension()` and `coerceToExtensionClass()`.
- Record the change in `NEWS.md`:

  ```markdown
  - Replaced obsolete session extension registration and dynamic S4 marker
    classes with mizer's simpler S3 extension mechanism. Stored parameter
    objects are now S3 objects with extension class vector
    `c("<pkg>", "MizerParams")`, removing the need for a load hook or an active
    binding.
  ```

#### 10. Document, test and check

```r
devtools::document()   # updates NAMESPACE and man/*.Rd
devtools::test()       # run test suite
devtools::check()      # verify full package build and vignettes
```

## Part 2 — other mizer changes that hit extension code

Work through these whatever the state of the S3 migration; none of them depends
on it.

### Working under both numerical schemes

Mizer's numerics are no longer a single scheme. The `second_order_w` slot selects
how the model is discretised in size — `flux` picks the advective reconstruction
and `bin_average` decides whether a size-dependent factor is integrated over its
bin or point-sampled at the left bin edge — and `project()`'s `method` argument
selects the discretisation in time, `"euler"` or one of the second-order
stepper `"second_order"`, an alias for `"tr_bdf2"`, or the superseded
`"predictor_corrector"`. Both default to the
first-order choice, which is why this is easy to miss: an extension that ignores
them looks correct, passes its own tests, and is then silently wrong for the
user who switches them on. Aim to work under all of them, and where you cannot,
say so out loud.

**In size.** The rules are in the `extend-mizer` skill, under "Respecting the
model's quadrature scheme": build derived quantities from the rate functions
rather than re-deriving them, use `encounter_kernel()` rather than
`pred_kernel()` inside the convolution, write summary integrals with
`sizeIntegral()`, and gate any array your own setter precomputes on
`params@second_order_w[["bin_average"]]`. Two consequences reach further than
that section:

- **A custom `Mort` or `Diffusion` rate must be the bin average when
  `bin_average` is on**, because the transport step differences it against the
  bin-averaged density — `getDiffusion()` and `getMort()` already return the bin
  average when the flag is on, because the arrays behind them are built that way
  at setup, so a function that wraps them inherits it and one that builds the
  rate from scratch does not. A custom `EGrowth` stays a point value at the bin
  boundary under both settings: growth velocities are never bin-averaged.
- **Do not reimplement the transport step or the egg balance.** They read the
  `flux` entry, so `getRequiredRDD()`, `findSteadyState()` and `project()` all
  agree about which reconstruction the model uses, and a hand-rolled copy in
  your package will not.

**In time.** The rules are in the `extend-mizer` skill, under "Components and
the time stepper": under either second-order stepper a component's
`dynamics_fun` — and a custom `resource_dynamics` function — is called **twice
per time step**, once as the predictor with the start-of-step rates and once as
the corrector with the midpoint rates, both times from the same `t` and the same
start-of-step state, and only the corrector's return value survives. What to
check in code that was written before those steppers existed:

- **Be a pure function of your arguments.** Compute the new state from the
  `n_other[[component]]` you were handed and return it. Anything that instead
  accumulates into a stored value — `<<-`, writing into `other_params(params)`,
  a counter, a log file — runs twice per step and is wrong by a factor of two
  under the second-order methods while looking fine under `"euler"`.
- **Use the `rates` argument; do not recompute rates from `n`.** Rebuilding them
  gives the corrector the start-of-step rates again, which silently drops the
  whole run back to first order rather than producing a visible error.
- **Integrate your own step to second order.** The corrector hands you midpoint
  values for `e_growth`, `mort`, `diffusion`, `rdd` and `resource_mort`, and
  start-of-step values for every other entry of the rate list. The interface
  asks for the new state, so your function is the integrator: `state + dt * f(state, rates)` stays first order under
  `"second_order"`. Solve the step exactly with the rates frozen, as
  `resource_semichemostat()` does, or take an RK2 step inside the function.
- **Tolerate being called out of order.** The rate functions are evaluated a
  second time at `t + dt` on the *predicted* state, so a rate function must not
  assume that time advances once per call or cache "the last `t` I saw".
- The existing rule about discontinuities still applies, and is the one thing a
  better stepper cannot fix: a rate that jumps as a function of abundance is not
  rescued by `"second_order"`.

**Test all four corners.** Running the package's own fixture under both
switches, rather than only on the default path, is covered in the
`upgrade-mizer-code` skill under "If your code is a package". An extension has
more to get wrong in time, so add the order check: halving `dt` should move the
answer by roughly a quarter rather than a half.

**If you cannot support one, warn.** The size scheme is a property of the object,
so the setup function is the place to say so:

```r
if (isTRUE(second_order_w(params)[["bin_average"]])) {
    signal_info(
        "second_order_w",
        paste("myExtension's detritus encounter is point-sampled at the left",
              "bin edge, so this model is not second order in w."),
        level = 1, severity = "warning", unhandled = "show"
    )
}
```

The time stepper is chosen per call rather than stored, so there is no equally
good hook. Name the limitation in the same setup-time report, and add a
`project.myExtension` method that inspects `method` before `NextMethod()` if you
want a per-call warning as well — remembering that it catches only the calls
that go through `project()`. The steady-state finders run the dynamics through
an internal stepper, so `tuneSteadyState(method = "second_order")` will not
reach it.

### Methods on the renamed rate-array accessors

Seventeen accessors that read a stored parameter or rate array back out of a
`MizerParams` object lost their `get` prefix: `getCatchability()` →
`catchability()`, `getExtMort()` → `ext_mort()`, and so on. For a *user* nothing
breaks — the `get` names are kept as aliases. For an *extension* they can break
silently, because the `get` name is no longer a generic of its own: it is now
bound to the same function object as the bare name (`getExtMort <- ext_mort`). A
method registered on the `get` name is therefore never called.

Rename the method to the bare name. These eleven bare names are generics, so a
method on them dispatches through both names:

`catchability`, `selectivity`, `pred_kernel`, `search_vol`, `intake_max`,
`metab`, `ext_mort`, `ext_encounter`, `maturity`, `repro_prop`,
`reproduction_level`

```r
# BEFORE                          # AFTER
getExtMort.myExtension <- ...     ext_mort.myExtension <- ...
```

The remaining six — `initial_effort`, `interaction_matrix`, `resource_dynamics`,
`resource_level`, `resource_rate`, `resource_capacity` — are plain functions with
no `UseMethod()` call, so there is nowhere to hang a method at all. If your
package relied on overriding one of them, that override is dead code: remove it
and open an issue asking for the bare name to be made a generic, rather than
leaving a method that looks live and is not.

A method exported with `@method getExtMort myExtension` still appears in
`NAMESPACE` and still passes `R CMD check`, so nothing will point this out for
you. Grep the package for method definitions whose name starts with `get`.

### Methods on the steady-state finders

`steady()` and `projectToSteady()` are superseded by `tuneSteadyState()` and
`findSteadyState()`, plus `projectUntilSettled()` for the trajectory. Unlike the
accessors above, the old names are **still generics with `MizerParams` methods**,
so `steady.myExtension` and `projectToSteady.myExtension` keep dispatching and
nothing breaks.

But the four generics are independent: `steady()`'s method calls an internal
worker, not `tuneSteadyState()`. A method on one name is not reached through the
other. So:

- to support the new names, **add** `tuneSteadyState.myExtension` and
  `findSteadyState.myExtension`; do not move the old methods;
- keep the old methods for as long as your users may call the old names;
- note the renamed arguments when writing the new methods: `tol` becomes
  `distance_tol`, `t_per` becomes `t_check`, and `return_sim` is gone — use
  `projectUntilSettled()` for the trajectory.

### Withdrawing a species parameter column

If your package adds a species parameter column and takes it away again when the
user switches the extension off, it almost certainly does that by reaching into
the slot, because until mizer 3.3.1 nothing else worked:

```r
params@species_params[["my_col"]] <- NULL   # old workaround — remove
```

The supported route is assignment, which also updates `given_species_params()`
and triggers the rebuild:

```r
species_params(params)$my_col <- NULL
```

A column mizer knows how to calculate comes back as a calculated value; a column
of your own is gone. This needs `mizer (>= 3.3.1)` in `DESCRIPTION`. See "Taking
a species parameter column away again" in the `create-extension-package` skill.

### The `gamma`/`f0` workaround can go

The defaults for `gamma` and `f0` are now measured on mizer's own reference
state: a rate function registered with `setRateFunction()`, the `ext_encounter`
array, and every `other_encounter()` contribution — a component's
`encounter_fun` included — are excluded from the measurement. Packages that add
encounter, such as a temperature scalar or a second resource, used to see
`gamma` and `f0` drift when the model was rebuilt, and worked around it by
pinning the current `gamma` as a given species parameter. Drop the workaround:

```r
given_species_params(params)$gamma <- NULL
```

The `gamma` mizer then calculates describes the species' baseline search volume,
which is what a dynamic modulation is meant to modulate. Where a species
genuinely has no predation encounter in the reference state —
`interaction_resource = 0`, or a kernel that does not overlap the resource —
supply `gamma` explicitly: mizer now raises an error naming the species instead
of deriving a nonsense value from the additive contribution.

### Tests that match messages or assert silence

Rephrased messages and calls that used to be quiet break tests in any package,
and the sweep for them is in the `upgrade-mizer-code` skill under "If your code
is a package". Two of them land specifically on extension code:

- `tuneSteadyState()` warns when a custom `resource_dynamics` has no matching
  `balance_<dynamics>()` function, so a test asserting silence on a package that
  supplies its own resource dynamics will fail. That is the point of the change:
  supply the missing `balance_<dynamics>()` function rather than silencing it.
- Anything the package emits **itself** should go through `signal_info()` rather
  than `message()` or `warning()`, so that it honours `info_level` and is
  collected with mizer's own reports — see the `create-extension-package` skill,
  "Telling the user what your package decided".

### Links to renamed articles and skills

The cheatsheet articles became guides, and the extension material was resplit.
Fix these wherever the package links them — `README`, vignettes, roxygen
`@seealso`:

| Old | New |
|---|---|
| `articles/creating-extension-packages.html` | `articles/guide-create-extension-package.html` |
| `articles/using-extension-packages.html` | `articles/guide-use-extension-packages.html` |
| `articles/extending-mizer.html` | `articles/guide-extend-mizer.html` |
| `vignette("cheatsheet-…")` | `vignette("guide-…")` |
| the `build-multispecies-model` skill | the `build-model` skill |

The website redirects the old addresses, but a `vignette()` call in your own code
or examples does not. Two of these rows are more than a rename: "Extending
mizer" and "Guide: Extending mizer" merged into a single guide covering the
mechanisms for changing mizer's dynamics, while everything that only matters once
you share an extension moved into the guide to creating an extension package.
Prose that paraphrases the old "Creating extension packages" article should be
re-read against the new one: its advice on marker classes was wrong by the end,
still telling you to write `setClass("myExtension", contains = "MizerParams")`.

## Finishing

Run the checklist for package authors in the `create-extension-package` skill
against the result, and its procedure for running mizer's test suite against
your subclass if the package is a dispatching extension. A migrated package is
the best case for that check, since the whole point of the migration is that
mizer's own behaviour still holds for your subclass.
