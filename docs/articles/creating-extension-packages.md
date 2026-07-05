# Creating a mizer extension package

## Overview

This vignette explains how to turn a mizer extension into a proper R
package. It is aimed at users who are already comfortable writing custom
rate functions or adding components with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
and who now want to share their extension with others or use it across
multiple projects.

There are two kinds of extension package:

- A **metadata-only** extension registers itself in `params@extensions`
  for record-keeping, but does not change how any mizer generic function
  behaves. [mizerStarvation](https://sizespectrum.org/mizerStarvation/)
  is an example: it adds starvation mortality via the `other_mort`
  pipeline, but it does not need to override any user-facing mizer
  functions.

- A **dispatching** extension additionally defines a new object type so
  that mizer’s generic functions (such as
  [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
  or
  [`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md))
  can be made to behave differently for models built with that
  extension. [mizerShelf](https://sizespectrum.org/mizerShelf/) is an
  example: it adds detritus and carrion components and overrides
  [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
  to include their biomasses in the result.

This vignette covers both kinds, working through concrete examples from
each package.

To get started quickly, clone or fork
[mizerExtensionTemplate](https://github.com/sizespectrum/mizerExtensionTemplate),
a minimal working package that illustrates all the mechanisms described
here with inline comments explaining each step.

## Why a package?

A plain R script works fine for a single project. An R package becomes
worthwhile when you want:

- **Reuse across projects**: install once,
  [`library()`](https://rdrr.io/r/base/library.html) everywhere.
- **Composability**: when two extension packages are loaded at the same
  time, mizer can arrange their contributions to generic functions in
  the right order automatically (see [Daisy-chaining with
  `NextMethod()`](#daisy-chaining-with-nextmethod) below).
- **Testing and documentation**: a package gives you a natural home for
  `testthat` tests, roxygen2 documentation and a pkgdown website.
- **Version tracking**: `params@extensions` records the version of each
  extension package used to build a model. If a collaborator opens your
  saved `MizerParams` object in a different session, mizer can warn them
  if the required package is missing or outdated.

## Metadata-only extensions: mizerStarvation

[mizerStarvation](https://github.com/sizespectrum/mizerStarvation) adds
starvation mortality — an extra per-capita mortality term that kicks in
when a fish’s energy balance is negative. It does this via the
`other_mort` argument in
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md):
supplying a function that
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md) calls
and adds to the mizer mortality at every time step. No mizer generic
function needs to be overridden.

### Registering in `.onLoad`

Every extension package, even a metadata-only one, should announce
itself to mizer when it is loaded. Place a `.onLoad` function in a file
such as `R/mizerMyExtension-package.R`:

``` r

.onLoad <- function(libname, pkgname) {
  mizer::registerExtension(pkgname, requirement = "owner/mizerMyExtension")
}
```

[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md)
adds the package name to the session’s extension chain. The
`requirement` string is a `pak` installation spec that mizer uses if it
needs to install the package automatically. For packages on CRAN you can
use a minimum version string such as `"1.2.0"` instead; for GitHub-only
packages use the `"owner/repo"` form
(e.g. `"sizespectrum/mizerStarvation"`). You can also specify a specific
branch or version of the package, using the same syntax that the
[pak](https://pak.r-lib.org/) package uses. The call is safe to repeat:
if the package is already registered (for example because the user
called
[`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
twice), it returns silently.

### Recording the extension in `params@extensions`

When your package creates or modifies a `MizerParams` object it should
copy the session’s registered extension chain into the `@extensions`
slot:

``` r

setStarvation <- function(params, starv_coef = 10) {
    # ... set up the rate function, species parameters, etc. ...
    params@extensions <- mizer::getRegisteredExtensions()
    params
}
```

[`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md)
returns the full chain that `.onLoad` hooks have built up. Storing this
in the object serves two purposes:

1.  **Reproducibility record.** When the object is saved with
    [`saveParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
    and later loaded with
    [`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md),
    mizer reads `@extensions` and warns if any required package is not
    installed or is too old.
2.  **Class coercion.** For dispatching extensions (see below),
    [`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
    uses `@extensions` to restore the correct S4 class automatically.

## Dispatching extensions: mizerShelf

[mizerShelf](https://github.com/sizespectrum/mizerShelf) adds two
dynamical components — detritus and carrion — to a mizer model. Beyond
just computing them, it also needs to change what certain user-facing
functions return: for example,
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
should include the detritus and carrion biomasses alongside the species
biomasses.

This section explains how to achieve that without breaking the standard
mizer behaviour, and without preventing other extension packages from
also modifying the same function.

### The problem with simply overwriting a function

Suppose you define a new
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
function in your package that adds your extra biomasses to the result.
That works as long as your package is the only one that modifies
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md).
But what if a second extension package also wants to add its own extra
components?

If both packages replace
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
whichever one was loaded last wins, and the other’s contribution is
silently lost. There is no way for the two packages to compose their
changes.

### How R dispatches to the right function

R has a built-in mechanism for exactly this situation. Every object has
a **class** attribute — a character string (or a vector of strings) that
labels what kind of thing it is. When you call a function like
`getBiomass(params)`, R looks at the class of `params` and searches for
a version of `getBiomass` whose name ends in `.` followed by that class
label, such as `getBiomass.mizerShelf`. If it finds one, it calls it. If
not, it tries the next class in the vector, and so on until it reaches
the base class and calls the default version.

This mechanism is called **S3 dispatch**, but you do not need to know
that term to use it. What matters practically is:

- Give your params objects a distinctive class label
  (e.g. `"mizerShelf"`).
- Define functions named `<genericname>.<classname>`
  (e.g. `getBiomass.mizerShelf`).
- Inside those functions, call
  [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) to pass
  control to the next class in the chain before or after your own
  modifications.

### Daisy-chaining with `NextMethod()`

[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) is what makes
multiple extension packages compose gracefully. Suppose the class of
`params` is `c("mizerFoo", "mizerShelf", "MizerParams")`, meaning
`params` is simultaneously of type `mizerFoo`, type `mizerShelf`, and
the base type `MizerParams`. Then calling `getBiomass(params)` proceeds
like this:

1.  R finds `getBiomass.mizerFoo` and calls it.
2.  `getBiomass.mizerFoo` calls
    [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).
3.  R finds `getBiomass.mizerShelf` and calls it.
4.  `getBiomass.mizerShelf` calls
    [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).
5.  R finds the standard mizer `getBiomass.MizerParams` and calls it.
6.  The standard result is returned to step 4, shelf biomasses are
    added, the result is returned to step 2, and `mizerFoo`’s biomasses
    are added on top.

Each extension in the chain sees and extends the result of all the
extensions below it. The chain grows automatically as packages are
loaded, so the user does not need to coordinate anything manually.

For this to work, **every method must call
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html)** so it does not
accidentally short-circuit the chain below it. The only exception is the
base mizer method at the bottom of the chain.

### Defining marker classes

To give your params objects a distinctive class label, you need to
define an **S4 marker class** — a formal class that extends
`MizerParams` but adds no new data. All extension-specific data lives in
`other_params(params)` or in component parameters; the class is just a
label.

Place these calls in a file such as `R/myextension-class.R`:

``` r

#' @export
setClass("mizerShelf", contains = "MizerParams")

#' @export
setClass("mizerShelfSim", contains = "MizerSim")
```

The params class name must match the name you pass to
[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md).
The sim class name must be the params class name with `"Sim"` appended.
`MizerSim` objects are coerced to the sim class automatically by
[`project()`](https://sizespectrum.org/mizer/reference/project.md) once
you record the extension chain in `params@extensions` (see below).

If your extension is designed to stack on top of another (say
`mizerBase`), inherit from that package’s class instead:

``` r

setClass("mizerOuter", contains = "mizerBase")
setClass("mizerOuterSim", contains = "mizerBaseSim")
```

### Registering in `.onLoad`

The `.onLoad` hook for a dispatching extension is the same as for a
metadata-only one:

``` r

.onLoad <- function(libname, pkgname) {
  mizer::registerExtension(pkgname, requirement = "owner/myExtensionPackage")
}
```

When
[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md)
is called, mizer prepends the extension to the session’s chain, giving
it the highest dispatch priority. Because R always loads dependency
packages before the package that depends on them, the dependent package
ends up outermost. For example, if `mizerOuter` depends on `mizerShelf`:

1.  R loads `mizerShelf`, its `.onLoad` fires → chain:
    `c(mizerShelf = "1.0.0")`
2.  R loads `mizerOuter`, its `.onLoad` fires → chain:
    `c(mizerOuter = "0.3.0", mizerShelf = "1.0.0")`

The class hierarchy `c("mizerOuter", "mizerShelf", "MizerParams")`
mirrors this chain, so dispatch proceeds in the right order
automatically.

### Writing methods that call `NextMethod()`

Here is `getBiomass.mizerShelf` from mizerShelf. It calls
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) first to get the
standard mizer result, then appends the detritus and carrion biomasses:

``` r

#' @method getBiomass mizerShelf
#' @export
getBiomass.mizerShelf <- function(object, ...) {
    params <- object
    b <- NextMethod()                       # standard species biomasses

    d_biomass <- sum(params@initial_n_pp *
                     params@dw_full * params@w_full)
    b <- c(b, Detritus = d_biomass)

    other <- params@initial_n_other
    scalar_other <- Filter(function(x) is.numeric(x) && length(x) == 1, other)
    if (length(scalar_other) > 0) b <- c(b, unlist(scalar_other))
    b
}
```

Because
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md)
calls
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
internally, this single override makes biomass plots include detritus
and carrion without any further changes.

Always register S3 methods in your package’s `NAMESPACE` file. The
roxygen2 `@method` tag does this for you automatically:

``` r

#' @method getBiomass mizerShelf
#' @export
getBiomass.mizerShelf <- function(object, ...) { ... }
```

### Replacing `setRateFunction()` with method dispatch

Users who write mizer extensions often start by replacing one of the
built-in rate functions with
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md):

``` r

myEncounter <- function(params, n, n_pp, n_other, t = 0, ...) {
    enc <- mizerEncounter(params, n = n, n_pp = n_pp, n_other = n_other, t = t, ...)
    enc + extraEncounter(params, n, n_pp, n_other, t, ...)
}
params <- setRateFunction(params, "Encounter", "myEncounter")
```

This works well for a single user’s workflow, but it is not composable:
if two extension packages both call
`setRateFunction(params, "Encounter", ...)`, whichever runs last
silently overwrites the other. When you turn your extension into a
package, **replace
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
calls with `project*` methods** for your marker class. These methods
participate in the daisy-chain via
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html), so two packages
can both modify the same rate without conflict.

#### The `project*` generics

Every standard mizer rate function has a corresponding S3 generic that
extension-aware projections call during
[`project()`](https://sizespectrum.org/mizer/reference/project.md).
Define a method for whichever rate your extension modifies:

| [`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md) key | S3 generic to override |
|----|----|
| `"Encounter"` | [`projectEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md) |
| `"FeedingLevel"` | [`projectFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md) |
| `"EReproAndGrowth"` | [`projectEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md) |
| `"ERepro"` | [`projectERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md) |
| `"EGrowth"` | [`projectEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md) |
| `"Diffusion"` | [`projectDiffusion()`](https://sizespectrum.org/mizer/reference/mizerDiffusion.md) |
| `"PredRate"` | [`projectPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md) |
| `"PredMort"` | [`projectPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md) |
| `"FMort"` | [`projectFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md) |
| `"Mort"` | [`projectMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md) |
| `"RDI"` | [`projectRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md) |
| `"RDD"` | [`projectRDD()`](https://sizespectrum.org/mizer/reference/projectRDD.md) |
| `"ResourceMort"` | [`projectResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md) |

#### Converting an existing custom rate function

Remove the
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
call from your constructor and define a method for your marker class
instead:

``` r

#' @method projectEncounter mizerMyExtension
#' @export
projectEncounter.mizerMyExtension <- function(params, n, n_pp, n_other,
                                              t = 0, ...) {
    enc <- NextMethod()
    enc + extraEncounter(params, n, n_pp, n_other, t, ...)
}
```

[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) replaces the
explicit call to
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md).
It passes control down the chain — first to any lower extension’s
`projectEncounter` method, and ultimately to
`projectEncounter.MizerParams`, which performs the standard mizer
calculation. Each extension in the chain adds its contribution on top of
the one below it, in load order.

Three rules:

1.  **Always call
    [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html)** — omitting
    it silently drops all contributions from lower extensions in the
    chain.
2.  **Keep the same argument signature** as the generic, and include
    `...` so extra arguments pass through.
3.  **Do not call
    [`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)**
    in your constructor for any rate that your package handles via a
    `project*` method. The two mechanisms are separate and should not be
    mixed for the same rate within an extension package.

#### How `setRateFunction()` and `project*` methods interact

A user who calls `setRateFunction(params, "Encounter", "myFn")` is
asking for their function to completely replace the encounter
calculation for that specific `params` object. Mizer honours this: when
`myFn` is set,
[`projectEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
is not called at all for the Encounter rate, so no extension package’s
`projectEncounter` method will run for that rate either.

This means that if a user applies
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
to a rate that your extension package modifies via
`projectEncounter.mizerMyExtension`, your method will be silently
bypassed for that object. It is worth documenting this limitation for
your users.

### Creating objects: the two commands

A constructor function that returns a `mizerShelf` object must end with
these two lines:

``` r

params@extensions <- mizer::getRegisteredExtensions()
params <- mizer::coerceToExtensionClass(params)
```

Here is how `newDetritusCarrionParams()` uses them in mizerShelf:

``` r

newDetritusCarrionParams <- function(species_params, ...) {
    params <- newMultispeciesParams(species_params, ...,
                                    resource_dynamics = "detritus_dynamics")
    # ... set up rate functions, components, colours ...

    params@extensions <- mizer::getRegisteredExtensions()
    params <- mizer::coerceToExtensionClass(params)
}
```

#### What `params@extensions <- getRegisteredExtensions()` does

[`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md)
returns the current session’s extension chain — everything that has been
registered via `.onLoad` hooks or with
[`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md)
once R started. Storing this in the params object is like stamping the
object with a bill of materials: it records exactly which extension
packages were active when the object was created.

Note that `@extensions` records the **full chain** that was active at
creation time, not just the outermost extension. If `mizerShelf` was the
only registered extension, `@extensions` will be
`c(mizerShelf = "1.0.0")`. If a further outer extension was also loaded,
both appear in the chain.

When the object is later loaded from disk with
[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md),
mizer reads `@extensions` to check that all the recorded extensions are
installed in the current session and warns the user if any are missing
or outdated.

#### What `coerceToExtensionClass(params)` does

At this point `params` is still a plain `MizerParams` object as far as R
is concerned. If you called `getBiomass(params)` now, R would call the
standard mizer `getBiomass.MizerParams` rather than
`getBiomass.mizerShelf`.

[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md)
reads `params@extensions`, finds the outermost extension in the object’s
own recorded chain that provides a dispatch class, and promotes the
object to that S4 class. In our example, `params` becomes
`c("mizerShelf", "MizerParams")` and R will dispatch to
`getBiomass.mizerShelf` automatically.
[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
calls this automatically when restoring an object from disk.

Note that coercion is driven by the object’s recorded chain, not by what
extensions happen to be loaded in the current session. An object created
with only `mizerShelf` registered will remain a `mizerShelf` object even
if `mizerOuter` is also loaded.

This step cannot be replaced with a simple
`class(params) <- "mizerShelf"`: because `MizerParams` is a formal S4
class, R enforces the class hierarchy strictly and the direct assignment
would fail.
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md)
uses the appropriate S4 machinery internally.

### What about `MizerSim` objects?

You do not need to call
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md)
yourself for `MizerSim` objects. When
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
creates its output it calls
[`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md),
which in turn calls
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md)
on the new sim object. Because the `params` slot inside the sim already
has `@extensions` set, mizer knows to promote the sim to `mizerShelfSim`
automatically.

This means that after:

``` r

sim <- project(NWMed_params, t_max = 3)
```

`sim` is already of class `mizerShelfSim`, and any method you have
defined for that class — such as `getBiomass.mizerShelfSim` — will be
dispatched automatically.

## Checklist for package authors

When building a dispatching extension package, verify the following:

Define both `setClass("<myExtension>", contains = "MizerParams")` and
`setClass("<myExtension>Sim", contains = "MizerSim")` with no new slots.

Call `mizer::registerExtension(pkgname, requirement = ...)` in
`.onLoad`.

End every constructor with
`params@extensions <- getRegisteredExtensions()` and
`coerceToExtensionClass(params)`.

Register every S3 method in `NAMESPACE` (via `@method` + `@export`).

Call [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) in every
method override.

For each rate the package modifies during projection, define a
`project*` method (e.g. `projectEncounter.mizerMyExtension`) rather than
calling
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md).
See [Replacing `setRateFunction()` with method
dispatch](#replacing-setratefunction-with-method-dispatch).

Store all extension-specific state in `other_params(params)` or in new
components created with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md),
never in new S4 slots.

Run mizer’s own test suite against an object of your subclass to check
that your overrides do not break core behaviour. See [Running mizer’s
test suite against your
subclass](#running-mizers-test-suite-against-your-subclass).

For metadata-only packages, only the second and third items apply (and
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md)
is not needed — just record
`params@extensions <- getRegisteredExtensions()`).

## Running mizer’s test suite against your subclass

A dispatching extension overrides mizer generics, so it is easy to
accidentally break some core mizer behaviour that you did not mean to
change. A powerful way to catch this is to run **mizer’s own test
suite** with the shared test fixture replaced by an object of your
subclass. If your overrides are faithful extensions of the base
behaviour, the great majority of mizer’s tests should still pass; the
failures that remain pinpoint exactly where your class diverges from a
plain `MizerParams` object (and are often legitimate — e.g. a test that
hard-codes the single-resource object structure).

Most of mizer’s tests build on a small shared fixture, `NS_params_small`
(and a simulation `NS_sim_small` derived from it), defined in
`tests/testthat/helper.R`. The idea is to turn that fixture into an
object of your subclass and then run the suite.

### Step 1: Get the mizer source

The test files are not part of the installed package, so clone the mizer
source and check out the tag or commit that matches your installed
version:

``` r

packageVersion("mizer")   # note this, then `git checkout` the matching tag
```

### Step 2: Convert the shared fixture to your subclass

Edit `tests/testthat/helper.R` and, immediately **after**
`NS_params_small` has been fully built (just before `NS_sim_small` is
created), insert code that replaces it with an object of your class.
Because
[`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
sources the helper inside mizer’s own namespace, attached packages are
not on the lookup path there, so:

- call [`library(yourPackage)`](https://rdrr.io/r/base/library.html)
  (this also fires your `.onLoad` and registers the extension), and
- qualify your own functions with `yourPackage::`.

For example, [mizerMR](https://github.com/sizespectrum/mizerMR) (which
adds multiple resources) converts the single resource into two resources
like this:

``` r

suppressMessages(library(mizerMR))
local({
    p1 <- NS_params_small
    rp <- data.frame(resource = c("Res A", "Res B"),
                     kappa = p1@resource_params$kappa / 2,
                     lambda = p1@resource_params$lambda, r_pp = 4,
                     n = p1@resource_params$n, w_min = min(p1@w_full),
                     w_max = p1@resource_params$w_pp_cutoff)
    strip <- function(m) { dimnames(m) <- NULL; m }
    ir <- p1@species_params$interaction_resource
    NS_params_small <<- suppressMessages(mizerMR::setMultipleResources(
        p1, resource_params = rp,
        resource_interaction = strip(cbind(ir, ir)),
        resource_capacity = strip(rbind(p1@cc_pp / 2, p1@cc_pp / 2)),
        resource_rate = strip(rbind(p1@rr_pp, p1@rr_pp)),
        initial_resource = strip(rbind(p1@initial_n_pp / 2, p1@initial_n_pp / 2))))
})
```

For your own extension, replace this block with a call to your
constructor or conversion function, so that `NS_params_small` becomes an
object of your S4 class. `NS_sim_small` is built from it on the next
line and will then be of your `...Sim` class automatically.

### Step 3: Skip the tests that manipulate the extension chain

A few of mizer’s test files deliberately test the extension-chain
machinery itself by calling
[`clearExtensionChain()`](https://sizespectrum.org/mizer/reference/clearExtensionChain.md)
and
[`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md).
The registered chain is **global session state**, and because testthat
runs all files in one R session, once those tests clear or replace the
chain the shared fixture is left orphaned: every later test that touches
it then fails in
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
with *“This object uses mizer extensions but no compatible extension
chain is registered”*. This is an artefact of the shared session, not a
problem with your extension, so exclude those files before running:

``` r

chain_tests <- c("test-extension-dispatch.R",
                 "test-registerExtensions.R",
                 "test-io.R")
file.remove(file.path("tests/testthat", chain_tests))
```

(You will restore them in step 5.)

### Step 4: Run the suite

Run the tests from the root of the mizer source tree. Use
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
rather than
[`testthat::test_dir()`](https://testthat.r-lib.org/reference/test_dir.html):
it calls `load_all()`, which is required because some mizer tests use
mizer’s internal (unexported) functions.

``` r

devtools::test()
```

### Step 5: Restore the source

The edits above are only for this experiment. Undo them with:

``` r

# from the mizer source root
system("git checkout tests/testthat")
```

### Interpreting the results

With the chain tests excluded, the only failures left are genuine
differences between your subclass and a plain `MizerParams` object.
Typical, expected ones include tests that assert the exact structure of
a single-resource object, or that index a resource array as if it were
one-dimensional. Any *other* failure — especially in a generic you
override — is worth investigating, as it usually means your method does
not faithfully extend the base behaviour.

## Upgrading objects across versions of your extension

As your extension evolves you may change where or how it stores its data
in a `MizerParams` object. Users who saved a model with an older version
of your package then need that model migrated to the new layout. mizer
upgrades the *core* slots itself (see
[`?validParams`](https://sizespectrum.org/mizer/reference/validParams.md)),
and it lets your extension hook into the same machinery so that your
migration runs automatically.

### Record a version stamp on the object

The `@extensions` slot records, for each extension in the chain, both
the requirement string used for dispatch and the version of the
extension package that the object conforms to. Write entries with
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md)
rather than assigning to the slot directly. Stamp the installed version
**only** when you create your component (or when you upgrade the
object); for ordinary modifications call
[`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md)
without a `version`, so the existing stamp is preserved:

``` r

# in your setup function, when the component is first created:
params <- recordExtension(params, "myExtension",
                          version = as.character(packageVersion("myExtension")))

# on later modifications, preserve the stamp:
params <- recordExtension(params, "myExtension")
```

Stamping only on create/upgrade is important: if an ordinary
modification re-stamped to the installed version, an object could claim
to be current and skip a migration it actually needs.

### Provide an `upgrade` method

mizer registers its core upgrades as methods of the S3 generic
[`utils::upgrade()`](https://rdrr.io/r/utils/upgrade.html). Register
your own method for your subclass and have it perform **only** your
migration. It must be idempotent, must **not** call
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html), and must
**not** touch the version stamp — mizer’s orchestrator re-stamps the
object after calling your method.

``` r

#' @exportS3Method utils::upgrade
upgrade.myExtension <- function(object, ...) {
    # Detect the old layout structurally and migrate it. Safe to run twice.
    if (!is.null(object@other_params$old_location)) {
        object@other_params$new_location <- object@other_params$old_location
        object@other_params$old_location <- NULL
    }
    object
}
```

### How it fires

[`needs_upgrading()`](https://sizespectrum.org/mizer/reference/needs_upgrading.md)
returns `TRUE` when the core mizer version is out of date **or** when
any extension’s recorded stamp is missing or older than the installed
package version. When a user runs the object through
[`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
(directly, or via
[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md),
[`project()`](https://sizespectrum.org/mizer/reference/project.md), a
setter, …) mizer’s orchestrator runs the core upgrade if needed and then
calls each out-of-date extension’s `upgrade` method in turn, re-stamping
each afterwards. A missing stamp counts as out of date, so objects
created before you adopted version tracking are migrated and stamped on
first use. Because your method is idempotent this is always safe.

Note that calling [`upgrade()`](https://rdrr.io/r/utils/upgrade.html) on
an object directly only dispatches to a single method and does not run
the full chain — `validParams(params)` (or
[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md))
is the entry point users should rely on.

## See also

- [mizerExtensionTemplate](https://github.com/sizespectrum/mizerExtensionTemplate)
  — a template package that puts all the mechanisms described here into
  a minimal, working package you can clone and adapt.
- [`vignette("extending-mizer", package = "mizer")`](https://sizespectrum.org/mizer/articles/extending-mizer.md)
  for the full menu of extension mechanisms (custom rate functions,
  components, subclassing, and more).
- [`vignette("using-extension-packages", package = "mizer")`](https://sizespectrum.org/mizer/articles/using-extension-packages.md)
  which explains the extension chaining from the package user’s
  perspective, and how to check or change which extensions are active in
  a session or in a model object.
- [`?registerExtension`](https://sizespectrum.org/mizer/reference/registerExtension.md)
- [`?getRegisteredExtensions`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md)
- [`?coerceToExtensionClass`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md)
- [`?recordExtension`](https://sizespectrum.org/mizer/reference/recordExtension.md)
