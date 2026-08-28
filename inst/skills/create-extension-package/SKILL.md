---
name: create-extension-package
description: >-
  Turn a working mizer extension into a shareable R package and maintain it.
  Use for packaging custom rates or components; registering an extension in
  .onLoad; choosing metadata-only or dispatching marker classes; chaining S3
  methods with NextMethod(); recording versions with recordExtension();
  coercing, bundling, testing and upgrading extension objects; or making user
  reports obey info_level. For implementing the underlying extension mechanisms
  use the extend-mizer skill; for using an existing package use the
  use-extension-packages skill.
---

# Creating a mizer extension package

## Overview

This is how to turn a mizer extension into a proper R package. It assumes you
are already comfortable writing custom rate functions or adding components with
`setComponent()` — the `extend-mizer` skill covers those — and now want to share
your extension with others or use it across several projects.

There are two kinds of extension package:

- A **metadata-only** extension registers itself in `params@extensions` for
  record-keeping, but does not change how any mizer generic function behaves.
  [mizerStarvation](https://sizespectrum.org/mizerStarvation/) is an
  example: it adds starvation mortality via the `other_mort` pipeline, but it
  does not need to override any user-facing mizer functions.

- A **dispatching** extension additionally defines a new object type so that
  mizer's generic functions (such as `getBiomass()` or `plotBiomass()`) can be
  made to behave differently for models built with that extension.
  [mizerShelf](https://sizespectrum.org/mizerShelf/) is an example: it
  adds detritus and carrion components and overrides `getBiomass()` to include
  their biomasses in the result.

Both kinds are covered here, working through concrete examples from each
package.

To get started quickly, clone or fork
[mizerExtensionTemplate](https://github.com/sizespectrum/mizerExtensionTemplate),
a minimal working package that illustrates all the mechanisms described here
with inline comments explaining each step.

## Why a package?

A plain R script works fine for a single project. An R package becomes
worthwhile when you want:

- **Reuse across projects**: install once, `library()` everywhere.
- **Composability**: when two extension packages are loaded at the same time,
  mizer can arrange their contributions to generic functions in the right order
  automatically (see [Daisy-chaining with `NextMethod()`] below).
- **Testing and documentation**: a package gives you a natural home for
  `testthat` tests, roxygen2 documentation and a pkgdown website.
- **Version tracking**: `params@extensions` records the version of each
  extension package used to build a model. If a collaborator opens your saved
  `MizerParams` object in a different session, mizer can warn them if the
  required package is missing or outdated.

## Metadata-only extensions: mizerStarvation

[mizerStarvation](https://github.com/sizespectrum/mizerStarvation) adds
starvation mortality — an extra per-capita mortality term that kicks in when a
fish's energy balance is negative. It does this through the `other_mort`
pipeline: `setStarvation()` registers the name of its rate function, and
`getMort()` calls every function registered there and adds the result to the
mortality rate at every time step. No mizer generic function needs to be
overridden.

That is the right mechanism for an extra mortality that depends on the state of
the model but carries no state of its own, and `other_mort()` is how you
register it:

```r
setStarvation <- function(params, starv_coef = 10) {
    # ... set up the species parameters the rate function needs ...
    other_mort(params)[["starvation"]] <- "starvMort"
    params
}
```

`other_encounter()` does the same for a contribution to the encounter rate. Two
neighbouring mechanisms are easy to confuse with these:

- If the contribution is a **fixed array** rather than a function of the model
  state, use `ext_mort()` or `ext_encounter()` instead. Those are cheaper and
  are what mizer's own external mortality and external encounter use.
- If the contribution comes from a **component with dynamics of its own** — a
  detritus or carrion pool that is itself updated each time step — create it
  with `setComponent()` and give the contribution as its `mort_fun` or
  `encounter_fun` argument. That entry then belongs to the component:
  `getComponent()` reports it and `removeComponent()` removes it, and
  `other_mort()` deliberately does not list it.

Do not assign into `params@other_mort` directly, even though older versions of
mizer left no alternative and mizerStarvation still does so. A bare slot
assignment skips the check that the name really is a function and that it does
not collide with a component's.

### Registering in `.onLoad`

Every extension package, even a metadata-only one, should announce itself to
mizer when it is loaded. Place a `.onLoad` function in a file such as `R/mizerMyExtension-package.R`:

```r
.onLoad <- function(libname, pkgname) {
  mizer::registerExtension(pkgname, requirement = "owner/mizerMyExtension")
}
```

`registerExtension()` adds the package name to the session's extension chain.
The `requirement` string is a `pak` installation spec that mizer uses if it
needs to install the package automatically. For packages on CRAN you can use a
minimum version string such as `"1.2.0"` instead; for GitHub-only packages use
the `"owner/repo"` form (e.g. `"sizespectrum/mizerStarvation"`). You can also specify a specific branch or version of the package, using the same syntax that the `{pak}` package uses. The call is safe to repeat: if the package is
already registered (for example because the user called
`devtools::load_all()` twice), it leaves the chain unchanged and rebuilds the
chain's dynamic marker classes if the reload removed any of them.

### Recording the extension in `params@extensions`

When your package creates or modifies a `MizerParams` object, record that your
extension has actually been applied. Stamp the installed package version when
the component is first created; on later modifications preserve the existing
stamp:

```r
setStarvation <- function(params, starv_coef = 10) {
    # ... set up the rate function, species parameters, etc. ...
    version <- if ("mizerStarvation" %in% names(params@extensions)) {
        NULL
    } else {
        as.character(utils::packageVersion("mizerStarvation"))
    }
    params <- mizer::recordExtension(params, "mizerStarvation",
                                     version = version)
    params
}
```

`recordExtension()` takes the installation requirement from the registered
chain, preserves all existing entries and version stamps, and adds this
extension in the correct dispatch position. Do not copy the entire active
registry into the object: another extension can be loaded without having been
applied to this particular model. As extension setup functions are applied,
their individual `recordExtension()` calls build the object's chain.

Storing this record serves two purposes:

1. **Reproducibility record.** When the object is saved with `saveParams()` and
   later loaded with `readParams()`, mizer reads `@extensions` and warns if any
   required package is not installed or is too old.
2. **Class coercion.** For dispatching extensions (see below), `readParams()`
   uses `@extensions` to restore the correct S4 class automatically.

### Taking a species parameter column away again

An extension that adds a species parameter column usually needs to remove it
when the user switches the extension off. Assign a table without the column and
mizer takes it out of both `species_params()` and `given_species_params()`:

```r
setStarvation <- function(params, starv_coef = 10) {
    if (all(starv_coef == 0)) {
        species_params(params)$starv_coef <- NULL   # withdraw the column
        return(params)
    }
    species_params(params)$starv_coef <- starv_coef
    # ... set up the rate function ...
}
```

Do **not** reach into `params@species_params` to do this. The rule is the same
for `given_species_params(params)$starv_coef <- NULL`: a column mizer knows how
to calculate comes back as a calculated value, and a column of your own — which
mizer has no way of recalculating — is gone. Before mizer 3.3.1 neither of these
worked, so a package that needs the behaviour has to require that version.

## Dispatching extensions: mizerShelf

[mizerShelf](https://github.com/sizespectrum/mizerShelf) adds two dynamical
components — detritus and carrion — to a mizer model. Beyond just computing
them, it also needs to change what certain user-facing functions return: for
example, `getBiomass()` should include the detritus and carrion biomasses
alongside the species biomasses.

This section explains how to achieve that without breaking the standard mizer
behaviour, and without preventing other extension packages from also modifying
the same function.

### The problem with simply overwriting a function

Suppose you define a new `getBiomass()` function in your package that adds your
extra biomasses to the result. That works as long as your package is the only
one that modifies `getBiomass()`. But what if a second extension package also
wants to add its own extra components?

If both packages replace `getBiomass()`, whichever one was loaded last wins, and
the other's contribution is silently lost. There is no way for the two packages
to compose their changes.

### How R dispatches to the right function

R has a built-in mechanism for exactly this situation. Every object has a
**class** attribute — a character string (or a vector of strings) that labels
what kind of thing it is. When you call a function like `getBiomass(params)`, R
looks at the class of `params` and searches for a version of `getBiomass` whose
name ends in `.` followed by that class label, such as `getBiomass.mizerShelf`.
If it finds one, it calls it. If not, it tries the next class in the vector,
and so on until it reaches the base class and calls the default version.

This mechanism is called **S3 dispatch**, but you do not need to know that term
to use it. What matters practically is:

- Give your params objects a distinctive class label (e.g. `"mizerShelf"`).
- Define functions named `<genericname>.<classname>` (e.g. `getBiomass.mizerShelf`).
- Inside those functions, call `NextMethod()` to pass control to the next class
  in the chain before or after your own modifications.

### Daisy-chaining with `NextMethod()`

`NextMethod()` is what makes multiple extension packages compose gracefully.
Suppose the class of `params` is `c("mizerFoo", "mizerShelf", "MizerParams")`,
meaning `params` is simultaneously of type `mizerFoo`, type `mizerShelf`, and
the base type `MizerParams`. Then calling `getBiomass(params)` proceeds like
this:

1. R finds `getBiomass.mizerFoo` and calls it.
2. `getBiomass.mizerFoo` calls `NextMethod()`.
3. R finds `getBiomass.mizerShelf` and calls it.
4. `getBiomass.mizerShelf` calls `NextMethod()`.
5. R finds the standard mizer `getBiomass.MizerParams` and calls it.
6. The standard result is returned to step 4, shelf biomasses are added, the
   result is returned to step 2, and `mizerFoo`'s biomasses are added on top.

Each extension in the chain sees and extends the result of all the extensions
below it. The chain grows automatically as packages are loaded, so the user
does not need to coordinate anything manually.

For this to work, **every method must call `NextMethod()`** so it does not
accidentally short-circuit the chain below it. The only exception is the base
mizer method at the bottom of the chain.

### Marker classes: do not define them yourself

Your params objects carry a distinctive class label, an **S4 marker class** that
extends `MizerParams` but adds no new data. All extension-specific data lives in
`other_params(params)` or in component parameters; the class is only a label.

**Do not create it with `setClass()`.** Mizer creates it for you when your
package is loaded. `registerExtension()` recognises your package as a
dispatching extension from the S3 methods it registers for its marker class —
`getBiomass.mizerShelf()` and the like — which it can see because registering an
S3 method does not require the class to exist yet. It then defines the class at
the right point in the hierarchy relative to whatever other extensions are
loaded.

That last part is why the static definition is worse than merely redundant.
`setClass("mizerShelf", contains = "MizerParams")` fixes your class as a direct
child of `MizerParams`, a sibling of every other extension, and seals it; mizer
cannot then re-parent it into the chain, so your package can no longer be
combined with another. Leaving the class to mizer is what lets two independently
developed extensions be loaded together in either order.

So `R/myextension-class.R` holds documentation and no code:

```r
#' mizerShelf marker classes
#'
#' S4 marker subclasses of MizerParams and MizerSim that enable S3 dispatch for
#' the methods defined in this package. They add no slots and are created by
#' mizer when the package is loaded, not by a setClass() call here.
#'
#' @name mizerShelf-class
#' @keywords internal
NULL
```

The class name is the name you pass to `registerExtension()`, and the sim class
is that name with `"Sim"` appended. `MizerSim` objects are coerced to the sim
class automatically by `project()` once you record the extension chain in
`params@extensions` (see below).

This also means there is nothing to do for an extension designed to stack on top
of another: mizer works out the inheritance from the chain, so you never name
the other extension's class.

### Registering in `.onLoad`

The `.onLoad` hook for a dispatching extension is the same as for a
metadata-only one:

```r
.onLoad <- function(libname, pkgname) {
  mizer::registerExtension(pkgname, requirement = "owner/myExtensionPackage")
}
```

When `registerExtension()` is called, mizer prepends the extension to the
session's chain, giving it the highest dispatch priority. Because R always
loads dependency packages before the package that depends on them, the dependent
package ends up outermost. For example, if `mizerOuter` depends on `mizerShelf`:

1. R loads `mizerShelf`, its `.onLoad` fires → chain: `c(mizerShelf = "1.0.0")`
2. R loads `mizerOuter`, its `.onLoad` fires → chain: `c(mizerOuter = "0.3.0", mizerShelf = "1.0.0")`

The class hierarchy `c("mizerOuter", "mizerShelf", "MizerParams")` mirrors
this chain, so dispatch proceeds in the right order automatically.

### Bundled data objects

If your package ships a ready-made `MizerParams` or `MizerSim` object in its
`data/` directory, you need one extra step in `.onLoad`. R's lazy-loading
mechanism delivers data objects exactly as they were serialised, without
calling `coerceToExtensionClass()`. If you do nothing, users get a plain
`MizerParams` with the wrong S4 class and method dispatch silently misfires.

The fix is to replace the lazy-loaded binding with an **active binding** that
calls `coerceToExtensionClass()` on every access. Because the coercion happens
at the moment the user touches the object, it always reflects the full
extension chain that is active at that point — including any other extension
packages the user loaded after yours.

Add this to your `.onLoad`, once per bundled params or sim object:

```r
.onLoad <- function(libname, pkgname) {
    mizer::registerExtension(pkgname, requirement = "owner/myExtensionPackage")
    if (exists("my_example_params", envir = asNamespace(pkgname), inherits = FALSE)) {
        ns  <- asNamespace(pkgname)
        raw <- get("my_example_params", envir = ns)   # capture the raw object once
        makeActiveBinding("my_example_params",
                          fun = function() mizer::coerceToExtensionClass(raw),
                          env = ns)
    }
}
```

After this, users can write `myPackage::my_example_params` (or simply
`my_example_params` after `library(myPackage)`) and always get a properly
classed object, with no extra steps on their part.

> **Why not just do `my_example_params <<- coerceToExtensionClass(...)` in
> `.onLoad`?**  That coercion runs when *your* package loads, at which point
> only your extension is in the chain. If the user then loads a second
> extension package, the chain grows, but the already-promoted object is never
> updated. `makeActiveBinding` avoids this by re-running the coercion on every
> access.



### Writing methods that call `NextMethod()`

Here is `getBiomass.mizerShelf` from mizerShelf. It calls `NextMethod()` first
to get the standard mizer result, then appends the detritus and carrion
biomasses:

```r
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

Because `plotBiomass()` calls `getBiomass()` internally, this single override
makes biomass plots include detritus and carrion without any further changes.

Always register S3 methods in your package's `NAMESPACE` file. The roxygen2
`@method` tag does this for you automatically:

```r
#' @method getBiomass mizerShelf
#' @export
getBiomass.mizerShelf <- function(object, ...) { ... }
```

### Replacing `setRateFunction()` with method dispatch

Users who write mizer extensions often start by replacing one of the built-in
rate functions with `setRateFunction()`:

```r
myEncounter <- function(params, n, n_pp, n_other, t = 0, ...) {
    enc <- mizerEncounter(params, n = n, n_pp = n_pp, n_other = n_other, t = t, ...)
    enc + extraEncounter(params, n, n_pp, n_other, t, ...)
}
params <- setRateFunction(params, "Encounter", "myEncounter")
```

This works well for a single user's workflow, but it is not composable: if two
extension packages both call `setRateFunction(params, "Encounter", ...)`,
whichever runs last silently overwrites the other. When you turn your extension
into a package, **replace `setRateFunction()` calls with `project*` methods**
for your marker class. These methods participate in the daisy-chain via
`NextMethod()`, so two packages can both modify the same rate without conflict.

#### The `project*` generics

Every standard mizer rate function has a corresponding S3 generic that
extension-aware projections call during `project()`. Define a method for whichever
rate your extension modifies:

| `setRateFunction()` key | S3 generic to override |
|------------------------|------------------------|
| `"Rates"` | `projectRates()` |
| `"Encounter"` | `projectEncounter()` |
| `"FeedingLevel"` | `projectFeedingLevel()` |
| `"EReproAndGrowth"` | `projectEReproAndGrowth()` |
| `"ERepro"` | `projectERepro()` |
| `"EGrowth"` | `projectEGrowth()` |
| `"Diffusion"` | `projectDiffusion()` |
| `"PredRate"` | `projectPredRate()` |
| `"PredMort"` | `projectPredMort()` |
| `"FMort"` | `projectFMort()` |
| `"Mort"` | `projectMort()` |
| `"RDI"` | `projectRDI()` |
| `"RDD"` | `projectRDD()` |
| `"ResourceMort"` | `projectResourceMort()` |

#### Converting an existing custom rate function

Remove the `setRateFunction()` call from your constructor and define a method
for your marker class instead:

```r
#' @method projectEncounter mizerMyExtension
#' @export
projectEncounter.mizerMyExtension <- function(params, n, n_pp, n_other,
                                              t = 0, ...) {
    enc <- NextMethod()
    enc + extraEncounter(params, n, n_pp, n_other, t, ...)
}
```

`NextMethod()` replaces the explicit call to `mizerEncounter()`. It passes
control down the chain — first to any lower extension's `projectEncounter`
method, and ultimately to `projectEncounter.MizerParams`, which performs the
standard mizer calculation. Each extension in the chain adds its contribution
on top of the one below it, in load order.

Three rules:

1. **Always call `NextMethod()`** — omitting it silently drops all contributions
   from lower extensions in the chain.
2. **Keep the same argument signature** as the generic, and include `...` so
   extra arguments pass through.
3. **Do not call `setRateFunction()`** in your constructor for any rate that
   your package handles via a `project*` method. The two mechanisms are
   separate and should not be mixed for the same rate within an extension
   package.

#### How `setRateFunction()` and `project*` methods interact

A user who calls `setRateFunction(params, "Encounter", "myFn")` is asking for
their function to completely replace the encounter calculation for that specific
`params` object. Mizer honours this: when `myFn` is set, `projectEncounter()`
is not called at all for the Encounter rate, so no extension package's
`projectEncounter` method will run for that rate either.

This means that if a user applies `setRateFunction()` to a rate that your
extension package modifies via `projectEncounter.mizerMyExtension`, your
method will be silently bypassed for that object. It is worth documenting this
limitation for your users.

#### Mizer's default parameters do not go through your method

Mizer calculates some species parameters by putting the model into a reference
state and measuring a rate in it. `get_gamma_default()` is the one to know
about: it gives each species a search volume coefficient of 1, puts a power-law
prey spectrum in front of it, and measures the available energy to find the
`gamma` that delivers the target feeding level `f0`. `get_f0_default()` does the
inverse.

Those measurements use `mizerEncounter()`, not `getEncounter()`, so they do not
dispatch through your `projectEncounter` method and are unaffected by a
`setRateFunction()` registration. They also exclude the `ext_encounter` array
and functions registered with `other_encounter()`, including a component's
`encounter_fun`. The `gamma` mizer calculates therefore describes the species'
*baseline* search volume on the reference resource — which is exactly what a
dynamic modulation such as a temperature scalar is meant to modulate. Do not
declare `gamma` as a given species parameter merely to protect it from your own
method or component; that is no longer necessary, and it costs the model the
ability to let `gamma` follow `f0`.

### Creating objects: the two commands

A constructor function that returns a `mizerShelf` object must end with these
two lines:

```r
params <- mizer::recordExtension(
    params, "mizerShelf",
    version = as.character(utils::packageVersion("mizerShelf")))
params <- mizer::coerceToExtensionClass(params)
```

Here is how `newDetritusCarrionParams()` uses them in mizerShelf:

```r
newDetritusCarrionParams <- function(species_params, ...) {
    params <- newMultispeciesParams(species_params, ...,
                                    resource_dynamics = "detritus_dynamics")
    # ... set up rate functions, components, colours ...

    params <- mizer::recordExtension(
        params, "mizerShelf",
        version = as.character(utils::packageVersion("mizerShelf")))
    params <- mizer::coerceToExtensionClass(params)
}
```

#### What `recordExtension()` does

`recordExtension()` adds the extension that created or modified this object,
taking its installation requirement from the session registry. Existing
entries keep their position and version stamps; a new entry is prepended so the
object's chain stays ordered outermost first. The result is a bill of materials
for the extensions actually applied to this model, not every extension package
that happened to be loaded when it was created.

The `version` stamp says which package version's object layout the new component
conforms to. A constructor supplies it because it has just created that
component. An ordinary modifier calls `recordExtension(params, "mizerShelf")`
without `version`, preserving the stamp already on the object. If several
extension setup functions are applied, their calls accumulate the full object
chain.

When the object is later loaded from disk with `readParams()`, mizer reads
`@extensions` to check that all the recorded extensions are installed in the
current session and warns the user if any are missing or outdated.

#### What `coerceToExtensionClass(params)` does

At this point `params` is still a plain `MizerParams` object as far as R is
concerned. If you called `getBiomass(params)` now, R would call the standard
mizer `getBiomass.MizerParams` rather than `getBiomass.mizerShelf`.

`coerceToExtensionClass()` reads `params@extensions`, finds the outermost
extension in the object's own recorded chain that provides a dispatch class,
and promotes the object to that S4 class. In our example, `params` becomes
`c("mizerShelf", "MizerParams")` and R will dispatch to
`getBiomass.mizerShelf` automatically. `readParams()` calls this automatically
when restoring an object from disk.

Note that coercion is driven by the object's recorded chain, not by what
extensions happen to be loaded in the current session. An object created with
only `mizerShelf` registered will remain a `mizerShelf` object even if
`mizerOuter` is also loaded.

This step cannot be replaced with a simple `class(params) <- "mizerShelf"`:
because `MizerParams` is a formal S4 class, R enforces the class hierarchy
strictly and the direct assignment would fail. `coerceToExtensionClass()` uses
the appropriate S4 machinery internally.

### What about `MizerSim` objects?

You do not need to call `coerceToExtensionClass()` yourself for `MizerSim`
objects. When `project()` creates its output it calls `MizerSim()`, which
in turn calls `coerceToExtensionClass()` on the new sim object. Because the
`params` slot inside the sim already has `@extensions` set, mizer knows to
promote the sim to `mizerShelfSim` automatically.

This means that after:

```r
sim <- project(NWMed_params, t_max = 3)
```

`sim` is already of class `mizerShelfSim`, and any method you have defined
for that class — such as `getBiomass.mizerShelfSim` — will be dispatched
automatically.

## Telling the user what your package decided

When your package makes a choice on the user's behalf — filling in a default,
adjusting an input, declining to carry out an instruction — report it through
mizer's own mechanism rather than with a plain `message()` or `warning()`. A
plain `message()` ignores `info_level`, is not collected with the other reports,
and is swallowed on the `species_params<-()` path. Four functions are exported
for this:

| Call | For |
|---|---|
| `signal_info()` | Any report about a choice you made or an input you adjusted |
| `with_info_level()` | Wrapping your entry point, so the reports raised inside it are collected and given together |
| `signal_not_recalculated()` | Your setter left a hand-set array alone |
| `default_info_level()` | The default for your own `info_level` argument |

Give every entry point that reports an `info_level` argument, forward it, and
wrap the body:

```r
newFooParams <- function(species_params, ...,
                         info_level = default_info_level()) {
    with_info_level(info_level = info_level, {
        params <- newMultispeciesParams(species_params,
                                        info_level = info_level, ...)
        if (is.null(species_params$foo_rate)) {
            signal_info("foo_rate", "No `foo_rate` provided, using 0.1.",
                        level = 1)
            params <- setComponent(params, ...)
        }
        params
    })
}
```

Handlers nest by themselves — an inner one steps aside and lets the outermost do
the reporting — so wrap without checking what your caller did.

Two arguments of `signal_info()` are worth getting right, because both wrong
choices are invisible at the call site:

- `severity` is not a judgement about how serious the report is, it is a fact
  about whether the report has to survive `suppressMessages()`. Use
  `"warning"` when the user asked for something that is not happening, and the
  default `"info"` when you are telling them about a choice you made.
  `species_params<-()` suppresses messages over its recalculation, so an
  `"info"` report raised under it never reaches the user.
- `unhandled` decides what happens when nothing is collecting, for instance when
  your setter is called directly. `"drop"` says nothing, which suits chatter that
  only makes sense as part of a report about a whole model. `"show"` reports it
  there and then, which is right when this may be all the user hears.

Do **not** hard-code the level in the call — `newMultispeciesParams(sp,
info_level = 0, ...)` looks like a way to keep your constructor quiet, but a
user who passes `info_level` themselves then gets `formal argument "info_level"
matched by multiple actual arguments`, because their value arrives through `...`
alongside yours.

Progress reports are the exception to all of this: they have to appear while the
work is happening, and collected reports are given at the end, so use a plain
`message()` for those.

## Checklist for package authors

When building a dispatching extension package, verify the following:

- [ ] Do **not** define the marker classes with `setClass()` — mizer creates
  `<myExtension>` and `<myExtension>Sim` when the package loads, and a statically
  defined class cannot be chained with other extensions. See
  [Marker classes: do not define them yourself].
- [ ] Call `mizer::registerExtension(pkgname, requirement = ...)` in
  `.onLoad`.
- [ ] For every `MizerParams` or `MizerSim` object bundled in `data/`, add a
  `makeActiveBinding()` call in `.onLoad` so the object is coerced to the
  correct extension class on access. See [Bundled data objects].

- [ ] End every constructor by calling `recordExtension()` with the installed
  package version, then `coerceToExtensionClass(params)`.
- [ ] Register every S3 method in `NAMESPACE` (via `@method` + `@export`).
- [ ] Call `NextMethod()` in every method override.
- [ ] For each rate the package modifies during projection, define a
  `project*` method (e.g. `projectEncounter.mizerMyExtension`) rather than
  calling `setRateFunction()`. See [Replacing `setRateFunction()` with method dispatch].
- [ ] Store all extension-specific state in `other_params(params)` or in
  new components created with `setComponent()`, never in new S4 slots.
- [ ] Add and remove any species parameter column of your own through
  `species_params(params)$my_col <- value` and
  `species_params(params)$my_col <- NULL`, never by writing into the
  `@species_params` slot. See
  [Taking a species parameter column away again].
- [ ] Register an extra mortality or encounter term that carries no state of its
  own with `other_mort()` or `other_encounter()`, never by assigning into
  `params@other_mort` or `params@other_encounter` directly. See
  [Metadata-only extensions: mizerStarvation].
- [ ] Report anything you tell the user with `signal_info()` inside a
  `with_info_level()`, never a bare `message()` or `warning()`, and give every
  entry point an `info_level = default_info_level()` argument that it forwards.
  See [Telling the user what your package decided].
- [ ] Run mizer's own test suite against an object of your subclass to check
  that your overrides do not break core behaviour. See
  [Running mizer's test suite against your subclass].

For metadata-only packages, only the `registerExtension()` call, the
bundled-data binding, the storage item and the reporting item apply, and
`coerceToExtensionClass()` is not needed. Their setup functions should still
call `recordExtension()`, stamping the package version when their component is
first created.

## Running mizer's test suite against your subclass

A dispatching extension overrides mizer generics, so it is easy to
accidentally break some core mizer behaviour that you did not mean to change. A
powerful way to catch this is to run **mizer's own test suite** with the shared
test fixture replaced by an object of your subclass. If your overrides are
faithful extensions of the base behaviour, the great majority of mizer's tests
should still pass; the failures that remain pinpoint exactly where your class
diverges from a plain `MizerParams` object (and are often legitimate — e.g. a
test that hard-codes the single-resource object structure).

Most of mizer's tests build on a small shared fixture, `NS_params_small` (and a
simulation `NS_sim_small` derived from it), defined in
`tests/testthat/helper.R`. The idea is to turn that fixture into an object of
your subclass and then run the suite.

### Step 1: Get the mizer source

The test files are not part of the installed package, so clone the mizer source
and check out the tag or commit that matches your installed version:

```r
packageVersion("mizer")   # note this, then `git checkout` the matching tag
```

### Step 2: Convert the shared fixture to your subclass

Edit `tests/testthat/helper.R` and, immediately **after** `NS_params_small` has
been fully built (just before `NS_sim_small` is created), insert code that
replaces it with an object of your class. Because `devtools::load_all()` sources
the helper inside mizer's own namespace, attached packages are not on the lookup
path there, so:

- call `library(yourPackage)` (this also fires your `.onLoad` and registers the
  extension), and
- qualify your own functions with `yourPackage::`.

For example, [mizerMR](https://github.com/sizespectrum/mizerMR) (which adds
multiple resources) converts the single resource into two resources like this:

```r
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

For your own extension, replace this block with a call to your constructor or
conversion function, so that `NS_params_small` becomes an object of your S4
class. `NS_sim_small` is built from it on the next line and will then be of your
`...Sim` class automatically.

### Step 3: Skip the tests that manipulate the extension chain {#skip-chain-tests}

A few of mizer's test files deliberately test the extension-chain machinery
itself by calling `clearExtensionChain()` and `registerExtensions()`. The
registered chain is **global session state**, and because testthat runs all
files in one R session, once those tests clear or replace the chain the shared
fixture is left orphaned: every later test that touches it then fails in
`validParams()` with *"This object uses mizer extensions but no compatible
extension chain is registered"*. This is an artefact of the shared session, not
a problem with your extension, so exclude those files before running:

```r
chain_tests <- c("test-registerExtensions.R",
                 "test-saveParams.R")
file.remove(file.path("tests/testthat", chain_tests))
```

(You will restore them in step 5.)

### Step 4: Run the suite

Run the tests from the root of the mizer source tree. Use `devtools::test()`
rather than `testthat::test_dir()`: it calls `load_all()`, which is required
because some mizer tests use mizer's internal (unexported) functions.

```r
devtools::test()
```

### Step 5: Restore the source

The edits above are only for this experiment. Undo them with:

```r
# from the mizer source root
system("git checkout tests/testthat")
```

### Interpreting the results

With the chain tests excluded, the only failures left are genuine differences
between your subclass and a plain `MizerParams` object. Typical, expected ones
include tests that assert the exact structure of a single-resource object, or
that index a resource array as if it were one-dimensional. Any *other* failure —
especially in a generic you override — is worth investigating, as it usually
means your method does not faithfully extend the base behaviour.

## Upgrading objects across versions of your extension

As your extension evolves you may change where or how it stores its data in a
`MizerParams` object. Users who saved a model with an older version of your
package then need that model migrated to the new layout. mizer upgrades the
*core* slots itself (see `?validParams`), and it lets your extension hook into
the same machinery so that your migration runs automatically.

### Record a version stamp on the object

The `@extensions` slot records, for each extension in the chain, both the
requirement string used for dispatch and the version of the extension package
that the object conforms to. Write entries with `recordExtension()` rather than
assigning to the slot directly. Stamp the installed version **only** when you
create your component (or when you upgrade the object); for ordinary
modifications call `recordExtension()` without a `version`, so the existing
stamp is preserved:

```r
# in your setup function, when the component is first created:
params <- recordExtension(params, "myExtension",
                          version = as.character(packageVersion("myExtension")))

# on later modifications, preserve the stamp:
params <- recordExtension(params, "myExtension")
```

Stamping only on create/upgrade is important: if an ordinary modification
re-stamped to the installed version, an object could claim to be current and
skip a migration it actually needs.

### Provide an `upgrade` method

mizer registers its core upgrades as methods of the S3 generic
`utils::upgrade()`. Register your own method for your subclass and have it
perform **only** your migration. It must be idempotent, must **not** call
`NextMethod()`, and must **not** touch the version stamp — mizer's orchestrator
re-stamps the object after calling your method.

```r
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

`needs_upgrading()` returns `TRUE` when the core mizer version is out of date
**or** when any extension's recorded stamp is missing or older than the
installed package version. When a user runs the object through `validParams()`
(directly, or via `readParams()`, `project()`, a setter, …) mizer's orchestrator
runs the core upgrade if needed and then calls each out-of-date extension's
`upgrade` method in turn, re-stamping each afterwards. A missing stamp counts as
out of date, so objects created before you adopted version tracking are migrated
and stamped on first use. Because your method is idempotent this is always safe.

Note that calling `upgrade()` on an object directly only dispatches to a single
method and does not run the full chain — `validParams(params)` (or
`readParams()`) is the entry point users should rely on.

## See also

- [mizerExtensionTemplate](https://github.com/sizespectrum/mizerExtensionTemplate)
  — a template package that puts all the mechanisms described here into a
  minimal, working package you can clone and adapt.
- the `extend-mizer` skill for the full menu of extension mechanisms — custom
  rate functions, external encounter and mortality, components, and subclassing.
- the [guide to using mizer extension
  packages](guide-use-extension-packages.html), which explains the extension
  chaining from the package user's perspective, and how to check or change
  which extensions are active in a session or in a model object.
- `?registerExtension`
- `?getRegisteredExtensions`
- `?coerceToExtensionClass`
- `?recordExtension`
