# Upgrading your extension package

This article is for the maintainer of an existing mizer extension
package — `mizerShelf`, `mizerMR`, `mizerStarvation`, `therMizer` or
your own — bringing it up to date with a newer mizer. It has two halves:

- **Part 1** is the one-off migration to mizer’s S3 architecture and
  simplified extension mechanism (mizer 3.4). Every package written
  against mizer 3.3 or earlier needs it.
- **Part 2** collects the other mizer changes that land specifically on
  extension code — methods that no longer dispatch, workarounds that can
  be dropped, tests that now fail. These are release-by-release and each
  one is independent of the S3 migration.

It covers only what is specific to *being an extension*.

------------------------------------------------------------------------

## Before you start: is this an extension?

Not every package that depends on mizer is a mizer extension, and only
an extension needs this article. Look for any of these:

- `setClass()` declaring a marker class, or `.onLoad()` calling
  `registerExtension()` / `registerExtensions()`
- methods on mizer’s own generics — `project.myClass`,
  `setParams.myClass`, `steadySingleSpecies.myClass`
- a rate function installed with
  [`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md),
  or a component registered with
  [`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)

A package with none of them — one that only *calls* mizer, however
extensively, like `mizerExperimental` — is ordinary user code in a
package wrapper. The [Upgrading your mizer code
article](https://sizespectrum.org/mizer/articles/upgrading.md) covers it
in full, including its “If your code is a package” section. Stop here
and use that.

**Everything in this article is in addition to that skill, not instead
of it.** Work through the [Upgrading your mizer code
article](https://sizespectrum.org/mizer/articles/upgrading.md) first —
its code-pattern index for the asymptomatic changes, its symptom index
for anything that broke, and its package section for the version floor,
the changelog, the docs and the test suite. Then come back here for what
only an extension has to do. For what the finished package should look
like, rather than how to get there, read the [guide to creating a mizer
extension
package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.md):
this article does not repeat its rules, and its checklist for package
authors is the acceptance test for a migrated package.

------------------------------------------------------------------------

## Part 1 — the S3 extension mechanism

### What changed in mizer

Starting with mizer 3.3.0.9000 / 3.4:

1.  **[`MizerParams`](https://sizespectrum.org/mizer/reference/MizerParams.md)
    and
    [`MizerSim`](https://sizespectrum.org/mizer/reference/MizerSim.md)
    are S3 classes.** They are stored as lists with class vectors
    `c("MizerParams")` and `c("MizerSim")`.
2.  **S3 extension classes:** an extension package prepends its class
    name to the S3 class vector — `c("mizerShelf", "MizerParams")` for
    params, and `c("mizerShelfSim", "MizerSim")` for simulations.
3.  **Session registration is obsolete:** `registerExtension()` and the
    dynamic S4 marker classes created at package load time have been
    removed.
4.  **No `.onLoad` hook or active bindings needed:** bundled datasets
    saved as S3 objects retain their class vectors and extension
    metadata under standard R lazy-loading.
5.  **Extension tracking via
    [`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md):**
    the extension requirement and package version are recorded in
    `params$extensions` (or `params@extensions`).
6.  **Class coercion via
    [`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md):**
    it reads `params$extensions` and sets the S3 class vector
    accordingly.
    [`project()`](https://sizespectrum.org/mizer/reference/project.md)
    automatically coerces its `MizerSim` output to the corresponding
    `…Sim` extension class.

### Step by step

#### 1. Clean up `R/<pkgname>-package.R`

Remove `.onLoad()` entirely (and delete `R/zzz.R` if it held only
`.onLoad()`), or the `registerExtension()` and
[`makeActiveBinding()`](https://rdrr.io/r/base/bindenv.html) calls
within it if it does anything else.

- **Remove `registerExtension()`**: session registration is no longer
  used.
- **Remove
  [`makeActiveBinding()`](https://rdrr.io/r/base/bindenv.html)**:
  datasets in `data/` now carry their S3 class directly.
- **Remove `@importFrom methods is`**: replace `is(x, "MizerParams")` in
  code with `inherits(x, "MizerParams")`.
- **Remove deleted files from `Collate:`**: if `DESCRIPTION` has a
  manual `Collate:` field, remove `zzz.R` (or whichever file held
  `.onLoad()`) from it.

``` r

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

An extension class is now an ordinary entry in the object’s S3 class
vector, so nothing has to be declared. Delete

``` r

setClass("myExtension", contains = "MizerParams")
setClass("myExtensionSim", contains = "MizerSim")
```

and, if the package carried the workaround for the marker classes being
wiped by `rm(list = ls())` or by the `cleanEx()` that `R CMD check` runs
before every example, delete that too:

``` r

setClass("myExtension", contains = "MizerParams", where = globalenv())
```

Both are actively harmful now: a sealed S4 class of the same name cannot
be re-parented, which is what prevented two independently developed
extensions from being chained in either load order.

Then check `NAMESPACE` and the roxygen comments for `@exportClass`,
`@importClassesFrom mizer` and `setValidity()` on those names, and
remove those as well.
[`devtools::document()`](https://rdrr.io/pkg/devtools/man/document.html)
will not drop a stale `exportClasses()` line while the `setClass()` is
still there.

#### 3. Update the class documentation

Rewrite the roxygen block that documented the extension class:

``` r

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

1.  call
    [`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md)
    with the package name, `version` and explicit
    `requirement = "<owner>/<repo>"` (now required since session
    registration no longer supplies it), and
2.  call `coerceToExtensionClass(params)`, for a dispatching extension.

``` r

params <- mizer::recordExtension(
    params, "myExtension",
    version = as.character(utils::packageVersion("myExtension")),
    requirement = "sizespectrum/myExtension"
)
params <- mizer::coerceToExtensionClass(params)
```

For a **metadata-only** extension such as `mizerStarvation`, call
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md)
and omit
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md).

#### 5. Migrate bundled data objects in `data/`

`MizerParams` and `MizerSim` objects in `data/` that were stored as S4
objects have to be upgraded to S3 and re-saved. Run this once, per
object:

``` r

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

``` r

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
- In `tests/testthat/`, replace `expect_s4_class()` with
  `expect_s3_class()` throughout — for `MizerParams` and `MizerSim` as
  well as for your own `myExtension` and `myExtensionSim`.
- Replace any assertions on `mizer::getRegisteredExtensions()` with
  assertions on `mizer::getMetadata(params)$extensions` (e.g., checking
  [`getMetadata(params)$extensions$myExtension$version`](https://sizespectrum.org/mizer/reference/setMetadata.md)
  and `$requirement`).

#### 7. Provide an `upgrade` method if layout has evolved

If your extension’s internal layout has changed across versions,
register an idempotent S3 method for
[`utils::upgrade()`](https://rdrr.io/r/utils/upgrade.html) so that older
saved objects are migrated automatically when validated. See “Upgrading
objects across versions of your extension” in the [guide to creating a
mizer extension
package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.md):

``` r

#' @exportS3Method utils::upgrade
upgrade.myExtension <- function(object, ...) {
    # Detect old layout and migrate idempotently
    object
}
```

#### 8. Update `DESCRIPTION`

- Remove `methods` from `Imports:` if it was there only for the S4
  checks.
- Require `mizer (>= 3.3.0.9000)`, or `mizer (>= 3.4.0)` once that is
  released.
- If `DESCRIPTION` has an explicit `Collate:` field, update it to remove
  deleted files (such as `zzz.R`) and include any new ones (such as
  `upgrade.R`).

#### 9. Update vignettes, narrative docs and `NEWS.md`

- Wherever a vignette describes S4 marker classes, `.onLoad()` or
  `registerExtension()`, replace it with S3 extension classes,
  [`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md)
  and
  [`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md).

- Record the change in `NEWS.md`:

  ``` markdown
  - Replaced obsolete session extension registration and dynamic S4 marker
    classes with mizer's simpler S3 extension mechanism. Stored parameter
    objects are now S3 objects with extension class vector
    `c("<pkg>", "MizerParams")`, removing the need for a load hook or an active
    binding.
  ```

#### 10. Document, test and check

``` r

devtools::document()   # updates NAMESPACE and man/*.Rd
devtools::test()       # run test suite
devtools::check()      # verify full package build and vignettes
```

------------------------------------------------------------------------

## Part 2 — other mizer changes that hit extension code

Work through these whatever the state of the S3 migration; none of them
depends on it.

### Working under both numerical schemes

Mizer’s numerics are no longer a single scheme. The
[`second_order_w`](https://sizespectrum.org/mizer/reference/second_order_w.md)
slot selects how the model is discretised in size — `flux` picks the
advective reconstruction and `bin_average` decides whether a
size-dependent factor is integrated over its bin or point-sampled at the
left bin edge — and
[`project()`](https://sizespectrum.org/mizer/reference/project.md)’s
`method` argument selects the discretisation in time, `"euler"` or one
of the second-order stepper `"second_order"`, an alias for `"tr_bdf2"`,
or the superseded `"predictor_corrector"`. Both default to the
first-order choice, which is why this is easy to miss: an extension that
ignores them looks correct, passes its own tests, and is then silently
wrong for the user who switches them on. That user need not do it on
purpose: a model created with defaults edition 3 or later is on the
second-order scheme in size from the start. Aim to work under all of
them, and where you cannot, say so out loud.

**In size.** The rules are in the [guide to extending
mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.md),
under “Respecting the model’s quadrature scheme”: build derived
quantities from the rate functions rather than re-deriving them, use
[`encounter_kernel()`](https://sizespectrum.org/mizer/reference/encounter_kernel.md)
rather than
[`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
inside the convolution, write summary integrals with
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md),
and gate any array your own setter precomputes on
`params@second_order_w[["bin_average"]]`. Two consequences reach further
than that section:

- **A custom `Mort` or `Diffusion` rate must be the bin average when
  `bin_average` is on**, because the transport step differences it
  against the bin-averaged density —
  [`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md)
  and [`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md)
  already return the bin average when the flag is on, because the arrays
  behind them are built that way at setup, so a function that wraps them
  inherits it and one that builds the rate from scratch does not. A
  custom `EGrowth` stays a point value at the bin boundary under both
  settings: growth velocities are never bin-averaged.
- **Do not reimplement the transport step or the egg balance.** They
  read the `flux` entry, so
  [`getRequiredRDD()`](https://sizespectrum.org/mizer/reference/getRequiredRDD.md),
  [`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
  and [`project()`](https://sizespectrum.org/mizer/reference/project.md)
  all agree about which reconstruction the model uses, and a hand-rolled
  copy in your package will not.

**In time.** The rules are in the [guide to extending
mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.md),
under “Components and the time stepper”: under either second-order
stepper a component’s `dynamics_fun` — and a custom
[`resource_dynamics`](https://sizespectrum.org/mizer/reference/setResource.md)
function — is called **twice per time step**, once as the predictor with
the start-of-step rates and once as the corrector with the midpoint
rates, both times from the same `t` and the same start-of-step state,
and only the corrector’s return value survives. What to check in code
that was written before those steppers existed:

- **Be a pure function of your arguments.** Compute the new state from
  the `n_other[[component]]` you were handed and return it. Anything
  that instead accumulates into a stored value — `<<-`, writing into
  [`other_params(params)`](https://sizespectrum.org/mizer/reference/setRateFunction.md),
  a counter, a log file — runs twice per step and is wrong by a factor
  of two under the second-order methods while looking fine under
  `"euler"`.
- **Use the `rates` argument; do not recompute rates from `n`.**
  Rebuilding them gives the corrector the start-of-step rates again,
  which silently drops the whole run back to first order rather than
  producing a visible error.
- **Integrate your own step to second order or supply an equation.** The
  corrector hands you midpoint values for all standard numeric rates. A
  custom step owns its integration: `state + dt * f(state, rates)`
  remains first order when `f` depends on the component’s start state.
  Solve the step exactly with the rates frozen, as
  [`resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.md)
  does, or use the optional prediction context described below. An
  equation-only component gets Heun integration from mizer.
- **Tolerate being called out of order.** The rate functions are
  evaluated a second time at `t + dt` on the *predicted* state, so a
  rate function must not assume that time advances once per call or
  cache “the last `t` I saw”.
- The existing rule about discontinuities still applies, and is the one
  thing a better stepper cannot fix: a rate that jumps as a function of
  abundance is not rescued by `"second_order"`.

**Test all four corners.** Running the package’s own fixture under both
switches, rather than only on the default path, is covered in the
[Upgrading your mizer code
article](https://sizespectrum.org/mizer/articles/upgrading.md) under “If
your code is a package”. An extension has more to get wrong in time, so
add the order check: halving `dt` should move the answer by roughly a
quarter rather than a half.

**If you cannot support one, warn.** The size scheme is a property of
the object, so the setup function is the place to say so:

``` r

if (isTRUE(second_order_w(params)[["bin_average"]])) {
    signal_info(
        "second_order_w",
        paste("myExtension's detritus encounter is point-sampled at the left",
              "bin edge, so this model is not second order in w."),
        level = 1, severity = "warning", unhandled = "show"
    )
}
```

The time stepper is chosen per call rather than stored, so there is no
equally good hook. Name the limitation in the same setup-time report,
and add a `project.myExtension` method that inspects `method` before
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) if you want a
per-call warning as well — remembering that it catches only the calls
that go through
[`project()`](https://sizespectrum.org/mizer/reference/project.md). The
steady-state finders run the dynamics through an internal stepper, so
[`tuneSteadyState(method = "second_order")`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
will not reach it.

### Your own FFT convolutions

Some extensions evaluate an encounter rate of their own — for a second
resource, a size-structured prey component, a diet split by prey group —
by copying the convolution from
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md):

``` r

avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) *
                                     mvfft(base::t(prey)),
                                 inverse = TRUE))) / length(params@w_full)
avail_energy <- avail_energy[, idx_sp, drop = FALSE]
```

That convolution is circular: the part of a small predator’s kernel that
reaches below the smallest size `min_w_pp` wraps around and picks up the
largest fish at the top of the grid (mizer issue \#618). From defaults
edition 3, mizer stores the Fourier-transformed kernels
`ft_pred_kernel_e`, `ft_pred_kernel_p` and `ft_pred_kernel_d` on a
zero-padded grid that is longer than
[`w_full`](https://sizespectrum.org/mizer/reference/w.md), so for a
model of that edition the code above fails with “non-conformable
arrays”. For a model of edition 1 or 2 the kernels keep their old length
and the code still runs, so the failure only shows up once a user works
with an edition-3 model. Replace it with
[`convolve_pred_kernel()`](https://sizespectrum.org/mizer/reference/convolve_pred_kernel.md),
which reads the length of the grid from the kernels and pads the input
to match:

``` r

avail_energy <- convolve_pred_kernel(params, prey)
```

It already returns only the predator sizes in `w`, so the `idx_sp`
subsetting goes. A vector `prey` is used for every predator species,
which is what a diet split by prey group needs: call it once per prey
group. For a predation rate, pass the predator densities on `w` with
`kernel = "predation"`; the result is on `w_full`. For a model of
edition 1 or 2 it gives the same circular convolution as the code it
replaces, so your extension keeps the results it had for those models.
It also uses the bin-integrated kernels when `bin_average` is on,
without any code of its own. This needs `mizer (>= 3.4.0.9002)` in
`DESCRIPTION`.

### The diet of a changed encounter rate

[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md) and
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md)
split the encounter rate between the food sources with the S3 generic
[`getEncounterByPrey()`](https://sizespectrum.org/mizer/reference/getEncounterByPrey.md),
and take the total from the model’s own encounter rate,
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md).
If your package changes the encounter rate with a
[`projectEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
method,
[`getEncounterByPrey()`](https://sizespectrum.org/mizer/reference/getEncounterByPrey.md)
does not know about the change, and the split no longer adds up to the
total. Mizer then scales the contributions of all food sources in
proportion and says so, with a message from
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`plotDiet()`](https://sizespectrum.org/mizer/reference/plotDiet.md) or
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md).
(Earlier versions ignored the change: the diet did not add up to what
the predator consumed, and the trophic levels could be far too high.)
The proportional split is exact only if your change affects all food
sources equally. Even a change of the search volume does not, as soon as
the model has an external encounter rate or a component with an
encounter function.

Give
[`getEncounterByPrey()`](https://sizespectrum.org/mizer/reference/getEncounterByPrey.md)
a method that makes the same change as your
[`projectEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
method and calls
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html):

``` r

#' @method projectEncounter myExtension
#' @export
projectEncounter.myExtension <- function(params, n, n_pp, n_other, t = 0, ...) {
    params <- withMySearchVolume(params, t)
    NextMethod()
}

#' @method getEncounterByPrey myExtension
#' @export
getEncounterByPrey.myExtension <- function(params, n, n_pp, n_other, t = 0,
                                           ...) {
    params <- withMySearchVolume(params, t)
    NextMethod()
}
```

If your change adds encounter of its own, a second resource for example,
add it to the result of
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) as a food source
of its own, named in the `prey` dimension, rather than spreading it over
the existing ones:

``` r

getEncounterByPrey.myExtension <- function(params, n, n_pp, n_other, t = 0,
                                           ...) {
    by_prey <- NextMethod()
    extra <- myResourceEncounter(params, n_other, t)  # species x size
    dn <- dimnames(by_prey)
    dn$prey <- c(dn$prey, "My resource")
    array(c(by_prey, extra), dim = dim(by_prey) + c(0, 0, 1), dimnames = dn)
}
```

[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md) then
gives a column for it, and
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
which only counts prey species and the resource in its numerator, treats
it like the external encounter rate.

[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md)
takes from your method only how much food each prey species and the
resource provide. It averages their trophic levels over prey sizes with
the predation kernel stored in the model, so a change to which prey
sizes are eaten reaches it only if the changed kernel is stored with
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md).
If your package changes the kernel during the projection instead, for
example with temperature, give
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md)
a method that stores the kernel in force at time `t` with
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
and calls [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).

A [`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md)
method of your own that exists only to account for your change to the
encounter rate is no longer needed: replace it with the
[`getEncounterByPrey()`](https://sizespectrum.org/mizer/reference/getEncounterByPrey.md)
method. If it calls
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) and then
computes the diet itself, that call now gives the message even though
its result is thrown away.
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md) also
gained a `t` argument, and
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md) on a
simulation passes the time of each saved step, so a method you keep must
accept `t`, at least through `...`, and pass it on. This needs
`mizer (>= 3.4.0.9003)` in `DESCRIPTION`.

### Methods on the renamed rate-array accessors

Seventeen accessors that read a stored parameter or rate array back out
of a `MizerParams` object lost their `get` prefix:
[`getCatchability()`](https://sizespectrum.org/mizer/reference/superseded_accessors.md)
→
[`catchability()`](https://sizespectrum.org/mizer/reference/setFishing.md),
[`getExtMort()`](https://sizespectrum.org/mizer/reference/superseded_accessors.md)
→
[`ext_mort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
and so on. For a *user* nothing breaks — the `get` names are kept as
aliases. For an *extension* they can break silently, because the `get`
name is no longer a generic of its own: it is now bound to the same
function object as the bare name (`getExtMort <- ext_mort`). A method
registered on the `get` name is therefore never called.

Rename the method to the bare name. These eleven bare names are
generics, so a method on them dispatches through both names:

`catchability`,
[`selectivity`](https://sizespectrum.org/mizer/reference/setFishing.md),
`pred_kernel`,
[`search_vol`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`intake_max`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`metab`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
`ext_mort`,
[`ext_encounter`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
[`maturity`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`repro_prop`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`reproduction_level`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)

``` r
# BEFORE                          # AFTER
getExtMort.myExtension <- ...     ext_mort.myExtension <- ...
```

The remaining six —
[`initial_effort`](https://sizespectrum.org/mizer/reference/initial_effort.md),
[`interaction_matrix`](https://sizespectrum.org/mizer/reference/setInteraction.md),
`resource_dynamics`,
[`resource_level`](https://sizespectrum.org/mizer/reference/setResource.md),
[`resource_rate`](https://sizespectrum.org/mizer/reference/setResource.md),
[`resource_capacity`](https://sizespectrum.org/mizer/reference/setResource.md)
— are plain functions with no
[`UseMethod()`](https://rdrr.io/r/base/UseMethod.html) call, so there is
nowhere to hang a method at all. If your package relied on overriding
one of them, that override is dead code: remove it and open an issue
asking for the bare name to be made a generic, rather than leaving a
method that looks live and is not.

A method exported with `@method getExtMort myExtension` still appears in
`NAMESPACE` and still passes `R CMD check`, so nothing will point this
out for you. Grep the package for method definitions whose name starts
with `get`.

### Methods on the steady-state finders

[`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.md)
and
[`projectToSteady()`](https://sizespectrum.org/mizer/reference/superseded_steady.md)
are superseded by
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
and
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md),
plus
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
for the trajectory. Unlike the accessors above, the old names are
**still generics with `MizerParams` methods**, so `steady.myExtension`
and `projectToSteady.myExtension` keep dispatching and nothing breaks.

But the four generics are independent:
[`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.md)’s
method calls an internal worker, not
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md).
A method on one name is not reached through the other. So:

- to support the new names, **add** `tuneSteadyState.myExtension` and
  `findSteadyState.myExtension`; do not move the old methods;
- keep the old methods for as long as your users may call the old names;
- note the renamed arguments when writing the new methods: `tol` becomes
  `distance_tol`, `t_per` becomes `t_check`, and `return_sim` is gone —
  use
  [`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
  for the trajectory.

### Withdrawing a species parameter column

If your package adds a species parameter column and takes it away again
when the user switches the extension off, it almost certainly does that
by reaching into the slot, because until mizer 3.3.1 nothing else
worked:

``` r

params@species_params[["my_col"]] <- NULL   # old workaround — remove
```

The supported route is assignment, which also updates
[`given_species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
and triggers the rebuild:

``` r

species_params(params)$my_col <- NULL
```

A column mizer knows how to calculate comes back as a calculated value;
a column of your own is gone. This needs `mizer (>= 3.3.1)` in
`DESCRIPTION`. See “Taking a species parameter column away again” in the
[guide to creating a mizer extension
package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.md).

### The `gamma`/`f0` workaround can go

The defaults for `gamma` and `f0` are now measured on mizer’s own
reference state: a rate function registered with
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md),
the `ext_encounter` array, and every
[`other_encounter()`](https://sizespectrum.org/mizer/reference/other_mort.md)
contribution — a component’s `encounter_fun` included — are excluded
from the measurement. Packages that add encounter, such as a temperature
scalar or a second resource, used to see `gamma` and `f0` drift when the
model was rebuilt, and worked around it by pinning the current `gamma`
as a given species parameter. Drop the workaround:

``` r

given_species_params(params)$gamma <- NULL
```

The `gamma` mizer then calculates describes the species’ baseline search
volume, which is what a dynamic modulation is meant to modulate. Where a
species genuinely has no predation encounter in the reference state —
`interaction_resource = 0`, or a kernel that does not overlap the
resource — supply `gamma` explicitly: mizer now raises an error naming
the species instead of deriving a nonsense value from the additive
contribution.

### Tests that match messages or assert silence

Rephrased messages and calls that used to be quiet break tests in any
package, and the sweep for them is in the [Upgrading your mizer code
article](https://sizespectrum.org/mizer/articles/upgrading.md) under “If
your code is a package”. Two of them land specifically on extension
code:

- [`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
  warns when a custom `resource_dynamics` has no matching
  `balance_<dynamics>()` function, so a test asserting silence on a
  package that supplies its own resource dynamics will fail. That is the
  point of the change: supply the missing `balance_<dynamics>()`
  function rather than silencing it.
- Anything the package emits **itself** should go through
  [`signal_info()`](https://sizespectrum.org/mizer/reference/signal_info.md)
  rather than [`message()`](https://rdrr.io/r/base/message.html) or
  [`warning()`](https://rdrr.io/r/base/warning.html), so that it honours
  `info_level` and is collected with mizer’s own reports — see the
  [guide to creating a mizer extension
  package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.md),
  “Telling the user what your package decided”.

### Links to renamed articles and skills

The cheatsheet articles became guides, and the extension material was
resplit. Fix these wherever the package links them — `README`,
vignettes, roxygen `@seealso`:

| Old | New |
|----|----|
| `articles/creating-extension-packages.html` | `articles/guide-create-extension-package.html` |
| `articles/using-extension-packages.html` | `articles/guide-use-extension-packages.html` |
| `articles/extending-mizer.html` | `articles/guide-extend-mizer.html` |
| `vignette("cheatsheet-…")` | `vignette("guide-…")` |
| the `build-multispecies-model` skill | the [guide to building a mizer model](https://sizespectrum.org/mizer/articles/guide-build-model.md) |

The website redirects the old addresses, but a
[`vignette()`](https://rdrr.io/r/utils/vignette.html) call in your own
code or examples does not. Two of these rows are more than a rename:
“Extending mizer” and “Guide: Extending mizer” merged into a single
guide covering the mechanisms for changing mizer’s dynamics, while
everything that only matters once you share an extension moved into the
guide to creating an extension package. Prose that paraphrases the old
“Creating extension packages” article should be re-read against the new
one: its advice on marker classes was wrong by the end, still telling
you to write `setClass("myExtension", contains = "MizerParams")`.

------------------------------------------------------------------------

## Stored component definitions

From mizer 3.4.0.9004,
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
accepts a reusable `definition` list and stores the assembled choices in
the model. This first preview supports the existing dynamics, encounter
and mortality callbacks and plotting styles. Continue to register
components through
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md);
direct edits to the legacy callback slots do not update the stored
definition. Those slots remain in use by the current runtime, including
temporary overrides during tuning.

Changing a definition object in your package affects subsequent
registrations, not models that have already used it. To migrate those
models, call
`setComponent(params, component, definition = new_definition)` from the
package’s upgrade method. Omit initial state and component parameters to
preserve each instance’s data. A replacement definition removes old
optional callbacks omitted from the new list; a call without
`definition` patches only supplied arguments.

The core upgrade creates snapshots for legacy components from their
recorded callbacks and plotting styles. It does not invent support for
future assessment, balancing or transformation capabilities. Keep using
the version-stamp and idempotent-upgrade workflow in the [guide to
creating a mizer extension
package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.md).

------------------------------------------------------------------------

## Component equations and prediction context

From mizer 3.4.0.9005, a component can register `rate_fun`,
`dynamics_fun`, or both. `rate_fun` returns a finite numeric derivative
with the state’s dimensions and names. Without a custom step, mizer uses
Euler or Heun; when both are registered, projection uses the custom step
and diagnostics can use the equation. Do not expect a fish TR-BDF2 step
to stabilise an explicit component equation. Keep a suitable custom step
for stiff or non-numeric components.

Custom component and resource steps can accept `predicted_state` (a list
of `n`, `n_pp`, `n_other` and `t`) and `endpoint_rates`. In a corrector
these are the predicted endpoint, while the ordinary state arguments and
`t` still refer to the start. In a predictor, Euler step or standalone
diagnostic these optional arguments are `NULL`. Strict legacy signatures
remain supported. All standard numeric fields of the corrector’s `rates`
list are now midpoint values; extra extension metadata is not
automatically averaged.

Declare a custom step’s `integration_order` only when its accuracy is
known. `sim@sim_params$component_integration` records the actual
per-component method, callback and order, with `NA` for an undeclared
custom order. Test temporal refinement with coupled states as well as a
component evolving on its own.

Callbacks may use exported `package::function` identifiers. The core
saved-model upgrade adds equation and method metadata without executing
callbacks. Existing models do not acquire equations merely because a
package’s definition changed; update their stored definitions explicitly
through the extension upgrade.

------------------------------------------------------------------------

## Registering convergence support

The assessment preview (3.4.0.9006) adds `convergence_fun` and `role` to
stored component definitions and
[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md)
for consumers and resource. Existing components upgrade as dynamic with
no invented convergence support. Projection continues to work, but an
operation requiring whole-model assessment can now report incomplete
coverage or stop before solving.

Supply a convergence callback returning named quantitative rows
(`measure`, `value`, `tolerance`, `units`). Values may be signed; mizer
tests their absolute values against the corresponding non-negative
tolerances. Different units can coexist. Provide meaningful measures for
non-biomass and non-numeric states too; do not wrap a maximum relative
cell rate merely to reproduce the old fallback. The callback receives
the full rate context and the current state and time. Strict callbacks
may accept only the context arguments they need.

For a custom resource, use
`setProcess(params, "resource", convergence_fun = "myResourceConvergence")`.
The core stores provider bindings, so a declaration cannot silently
follow a different dynamics implementation.
[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md)
is an explicit adoption of the current provider. When changing a
component’s equation or custom step, explicitly select the convergence
callback and role again if they still apply. A complete replacement
definition must include all the optional capabilities it is to retain.

Use `role = "prescribed"` only when the state block is an input of the
model. Freezing a dynamic state for one analysis instead uses explicit
conditional scope. Temporary frozen callbacks do not alter the stored
structural role. The built-in constant resource is recognised as a
prescribed input. The stability preview adds partial masks as described
below.

Update saved models through the package’s upgrade method to install
these declarations. Test passing, failing and missing coverage; signed
and zero-state measures; callback failures; provider replacement; and
both quadrature schemes for any size integrals. The current Newton
solver requires unsupported dynamic components to be explicitly fixed. A
convergence callback alone does not add their states to that solver.

------------------------------------------------------------------------

## Declaring structural inputs and temporal support

The stability preview (3.4.0.9010) stores independent
`prescribed_mask_fun`, `autonomous` and resource `rate_fun`
capabilities. Changing one capability does not renew the provider
bindings of the others. A provider replacement therefore needs explicit
adoption of each declaration that still applies.

A prescribed-mask callback receives the model, with optional `owner`,
`component` or `process` identifiers. It returns a logical mask with
exactly the state’s dimensions and names and no missing entries. `TRUE`
means the coordinate is a structural input, not an unknown in
equilibrium or stability analysis. The callback describes model
structure: do not infer a mask from a currently zero abundance or
derivative. Partial masks require a numeric state; `role = "prescribed"`
remains the shorthand for an entire input block.

Declare `autonomous = TRUE` only when the effective providers used by
that state owner have no explicit time dependence. This is separate from
`stationary` tuning support: `sin(t) * x` has a stationary solution at
zero but a time-dependent Jacobian there. Custom reproduction
contributes to the consumer equations and therefore needs the consumers’
temporal declaration. An explicit conditional exclusion removes the
excluded owner’s unused dynamics, but its encounter or mortality
contributions can still require temporal support.

For continuous resource stability, register `rate_fun` with
`setProcess(params, "resource", ...)`. It returns the actual resource
derivative with the same shape as
[`initialNResource(params)`](https://sizespectrum.org/mizer/reference/initialNResource-set.md)
and receives the complete effective rates, current state and time. Keep
a custom resource step as well when needed for projection. The built-in
semichemostat, logistic and constant resource laws have exact equation
adapters; mizer no longer estimates an unknown resource equation by
taking a very short custom step.

This preview’s Jacobians cover consumer and resource coordinates only.
Component equations and convergence callbacks do not automatically grant
Jacobian support. Whole-model requests with dynamic components report
missing coverage; callers must explicitly choose conditional scope to
hold them fixed. Test time-dependent forcing, partial prescribed masks,
zero-state invasion modes, stale bindings and unsupported drivers, as
well as the supported path.

------------------------------------------------------------------------

## Registering tuning and reproduction support

The tuning preview (3.4.0.9007) adds `stationary`, `equation_rates_fun`
and named `tuning` strategies to component and process definitions.
Upgrade legacy models by explicitly installing your package’s
declarations; the core upgrade does not infer them from a projecting or
constant-looking callback. The complete contract and field examples are
in the [guide to extending
mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.md)’s
“Declaring tuning strategies” section.

Each strategy partitions its owned states into held and evolving states,
declares permitted parameter paths and bounds, and supplies any
preserved quantities and author-written balancer. Declare a no-action
strategy explicitly when no parameters need adjustment. Coupled
responsibilities or shared adjustable parameters require a joint
strategy registered with
[`setTuningStrategy()`](https://sizespectrum.org/mizer/reference/setTuningStrategy.md).

Balancing must preserve the target state and every effective equation
rate. Report source and loss terms independently, rather than hiding
their changes inside a cancelling residual. Mizer constructs targets
before invoking the balancer and verifies the returned model
independently. It rejects undeclared parameter changes and has no
generic optimisation fallback or repair solve. Test infeasible targets
and parameter bounds, unselected species, zero states, and a
deliberately invalid balancer that preserves only the net residual.

Reproduction is a process without a separate state. Register its support
with `setProcess(params, "reproduction", ...)`, including a law-specific
`level_fun` when meaningful. The built-in law adapters apply only to
their actual effective providers; a named replacement or S3 wrapper
needs explicit support. In particular, a nonlinear or coupled RDI cannot
inherit an assumption that RDI is proportional to `erepro`. A custom
full-rates provider must expose the effective recruitment it supplies to
projection.

[`balanceRDD()`](https://sizespectrum.org/mizer/reference/balanceRDD.md)
wraps the method chain in one verification transaction. Methods can call
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) and then apply
their own declared adjustments, but the final result must reproduce the
effective target including additions such as immigration. Tuning pins
that final assembled recruitment after the complete rate pipeline;
intermediate RDI may change during balancing. Incidental calibration now
calls
[`balanceRDD()`](https://sizespectrum.org/mizer/reference/balanceRDD.md)
and retains the current reproduction law; reserve
[`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
for an intentional law switch.

------------------------------------------------------------------------

## Registering structural and scaling transformations

The transformation preview (3.4.0.9008) adds a `transformations` map to
component and process definitions. Legacy extensions do not acquire
support by having arrays whose dimensions happen to match species or
size grids. Declare an updating hook or explicit no-action support for
each operation the package supports; missing support stops the operation
before changes.

Hooks receive a candidate with the core structural bookkeeping
completed, plus the original model, operation details and ownership
context. They run before dependent rate construction, initialisation and
steady-state calculations. Move repairs such as a species-by-size
spatial-overlap array into these hooks, so
[`addSpecies()`](https://sizespectrum.org/mizer/reference/addSpecies.md)
does not evaluate rates while that array still has old dimensions. Do
not keep the same repair both in a hook and after
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) returns.

The public transformation functions retain their S3 generics. A
delegating method should let the coordinated hooks perform its data
updates, rather than doing the same update before or after
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html). A complete
overriding method remains responsible for implementing the operation’s
contract itself.

Use explicit species and grid mappings from `details`; choose how new
cells are initialised. For scaling, distinguish extensive stocks from
intensive quantities and non-numeric state.
[`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.md)
no longer multiplies arbitrary component state, and
[`scaleRates()`](https://sizespectrum.org/mizer/reference/scaleRates.md)
needs the component’s time constants and other rate parameters to
transform consistently too. Update species coefficients and their
recorded supplied values together where appropriate, and preserve the
intended freeze status of manually supplied arrays.

Declare permitted update paths. A component’s own state and parameter
subtree are its default paths; other extension-owned fields need
explicit paths. Shared data require joint ownership through
[`setTransformation()`](https://sizespectrum.org/mizer/reference/setTransformation.md).
Callback errors and undeclared changes are errors; an explicit no-action
declaration must really mean the component’s existing data remain valid
after the operation.

Support remains bound to the providers it was registered against.
Updating an assessment, a biomass measurement or instance data does not
re-adopt an old transformation or balancing declaration after a provider
change. Explicitly register the affected capability family, or replace
the complete definition, after checking that it supports the new
provider. A bare `setProcess(params, process)` deliberately re-adopts
the stored process contract as a whole. Reproduction-level measurement
has its own binding and cannot re-enable a stale balancing strategy.

For
[`removeSpecies(compensate = TRUE)`](https://sizespectrum.org/mizer/reference/removeSpecies.md),
the core adds each effective external-rate shortfall once and checks the
remaining species’ growth, diffusion, total mortality and recruitment
after all hooks. If your rate callback changes how those arrays enter
the equations, supply the corresponding correction in a removal hook
with explicit permission to update them. There is no numerical probe
that guesses your callback’s response. Custom resource dynamics and
other components own their compensation through the same hooks.

See “Transforming extension data with the model” in the [guide to
extending
mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.md)
for field names, callback context and a selective scaling example. Test
add/remove/rename/ reorder/grid changes and both scaling operations.
Include a dependent callback that checks dimensions when it is first
invoked, missing support, hook failure, unchanged input on failure, and
a saved-model round trip. Keep the public `species_params<-()` rule that
table rows alone cannot reorder the model; use
[`reorderSpecies()`](https://sizespectrum.org/mizer/reference/reorderSpecies.md)
for a coordinated reorder.

------------------------------------------------------------------------

## Registering biomass and mortality diagnostics

The measurement preview (3.4.0.9009) adds `biomass_fun` and
`biomass_filters` to component definitions. Register measurements
independently of convergence: even a component with a known biomass
needs its own stationarity assessment. The callback receives the current
state and time and returns named numeric biomasses or the richer
data-frame form documented in the [guide to extending
mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.md).
Declare supported filters explicitly and return unavailable rows with
reasons when a requested measurement cannot be made. Keep identifiers
and plotting metadata stable across simulation times, including when a
value is unavailable.

Move component discovery and biomass calculation from custom plotting
methods into the getter capability.
[`getBiomass(include = ...)`](https://sizespectrum.org/mizer/reference/getBiomass.md)
supplies labels, colours and coverage to
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md); plotting
should consume this result. The default stays species only. Do not make
a plot discover extra components behind the getter’s back, or sum
incomplete measurements with `na.rm = TRUE` and call that a total. Use
physical-source identifiers to prevent double counting a resource
exposed through more than one registration path.

For mortality plots, consume
[`getMortByCause()`](https://sizespectrum.org/mizer/reference/getMortByCause.md)
and its cause metadata. Registered additive mortality contributions
appear automatically with the built-in assembler. A replacement `Rates`
or `Mort` provider supplies a `mortality_fun` on the consumers process.
The callback’s complete decomposition must reproduce the effective total
at every state; an incomplete report keeps its missing-support reasons.
Adopting this diagnostic does not re-adopt stale tuning, convergence or
transformation declarations. Remove copied mortality plot logic only
after checking the extension’s intended causes and presentation against
the shared getter.

------------------------------------------------------------------------

## Reusing built-in capability definitions

The migration preview (3.4.0.9011) adds
[`mizerMortalityBreakdown()`](https://sizespectrum.org/mizer/reference/mizerMortalityBreakdown.md)
for an author callback to reuse the built-in terms before appending
extension causes. The same full effective rates are used, and
[`getMortByCause()`](https://sizespectrum.org/mizer/reference/getMortByCause.md)
checks the combined sum afterwards. Prefer this to a copied
decomposition or a separate call to
[`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md).
Predator-specific refinement is optional and must agree with the
effective aggregate; incomplete attribution never licenses inventing a
cause or renormalising known terms to the total.

For a known built-in reproduction law whose semantics your extension
preserves,
[`mizerReproductionDefinition()`](https://sizespectrum.org/mizer/reference/mizerReproductionDefinition.md)
supplies the corresponding strategies and level definition. This is an
explicit reuse tool, not automatic approval of a custom RDI, RDD or
complete-rate wrapper. Supply transformation declarations separately and
test preserved effective recruitment after balancing.

------------------------------------------------------------------------

## Registering capabilities during an upgrade

Extension upgrade methods can now register capabilities through
validating public setters without restarting their own migration. Remove
workarounds that stamp the extension as current before the upgrade has
finished. Mizer retains the old stamp during the callback, processes the
chain from inside out and records success afterwards.

------------------------------------------------------------------------

## Finishing

Run the checklist for package authors in the [guide to creating a mizer
extension
package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.md)
against the result, and its procedure for running mizer’s test suite
against your subclass if the package is a dispatching extension. A
migrated package is the best case for that check, since the whole point
of the migration is that mizer’s own behaviour still holds for your
subclass.
