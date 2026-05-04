# Extension Chain Implementation Plan

## Status

Replacement implementation plan for moving mizer extensions from
`rates_funcs`-style single-function overrides to a composable extension chain.

The design keeps `MizerParams` and `MizerSim` as the core S4 containers, but
uses dynamically registered S4 subclasses to expose an S3 class order. Rate
calculation hooks then become S3 generics, and extension packages compose by
calling `NextMethod()`.

## Key Decisions

- There is one maximal extension chain per R session.
- The extension vector order is the S3 dispatch order, from outermost to
  innermost.
- A registered session can also handle suffix subchains of the maximal chain.
- Base mizer objects with no extensions remain valid in an extension session.
- `MizerSim` does not get an `@extensions` slot. Its chain is always read from
  `sim@params@extensions`.
- Extension subclasses are marker classes for dispatch only. They must not add
  extra slots.
- Saved objects are serialised as base `MizerParams` or `MizerSim` objects so
  they can be read before dynamic subclass definitions exist.
- `rates_funcs` remains during the transition, but the new S3 dispatch path is
  the architecture extension packages should use.

## Extension Vector Semantics

An extension chain is a named character vector:

```r
c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")
```

The names are logical extension identifiers and params-side S4 class names. The
values are package requirements or installation specifications.

The order is left-to-right dispatch order:

```r
c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")
```

means:

1. call `projectEncounter.mizerExtB()`
2. then `NextMethod()` calls `projectEncounter.mizerExtA()`
3. then `NextMethod()` calls `projectEncounter.MizerParams()`

To achieve this, class definitions are built from right to left:

```text
mizerExtB contains mizerExtA contains MizerParams
```

so the in-memory object class is `mizerExtB`, and the S3 class vector seen by
`UseMethod()` is effectively:

```r
c("mizerExtB", "mizerExtA", "MizerParams")
```

## Subchains

The active session chain is the maximal chain. A valid object chain is either:

- identical to the active session chain,
- a suffix of the active session chain,
- or empty.

For example, if the session chain is:

```r
c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")
```

then these object chains are valid in the same session:

```r
c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")
c(mizerExtA = "0.4.1")
character()
```

but these are not valid:

```r
c(mizerExtB = "1.2.0")
c(otherExt = "1.0.0", mizerExtA = "0.4.1")
```

This supports the intended use case where one model uses an extra outer
extension and another model in the same session uses only the common inner
extension, or no extensions.

`registerExtensions()` may also grow the session chain by prepending outer
extensions. For example:

```r
registerExtensions(c(mizerExtA = "0.4.1"))
registerExtensions(c(mizerExtB = "1.2.0", mizerExtA = "0.4.1"))
```

is allowed because the first chain is a suffix of the second. Appending,
reordering, removing inner extensions, or changing version/source values is an
error.

## Session State

Add a package-private environment:

```r
.mizerSession <- new.env(parent = emptyenv())
```

with these fields:

- `.mizerSession$extensions`: the maximal active chain.

The empty session chain is represented by `character()`.
All params-side and sim-side class names are derived from this vector when
needed. No separate cached class lists are required.

## Class Naming

Params-side class names are the extension identifiers:

```r
mizerExtA
mizerExtB
```

Sim-side class names must be distinct from params-side class names because S4
class names are global. Use a deterministic helper:

```r
simExtensionClass <- function(extension) {
    paste0(extension, "Sim")
}
```

For `mizerExtA`, the sim class is `mizerExtASim`.

Extension classes must not add slots. They exist only to provide S3 dispatch
classes. Extension-specific persistent state should live in existing mizer
slots, such as `other_params`, `other_dynamics`, `other_encounter`,
`other_mort`, metadata, or other agreed extension storage. Adding slots would
make coercion from base objects, suffix subchains, save/load, and old-object
upgrade semantics much more complicated.

## Public And Internal API

### `registerExtensions()`

Suggested API:

```r
registerExtensions <- function(extensions, install = FALSE) {
}
```

Responsibilities:

- validate the extension vector.
- check required extension namespaces are available.
- optionally install missing extension packages if `install = TRUE`.
- load extension namespaces so S3 methods are registered.
- define params and sim S4 subclass chains.
- update the maximal session chain.
- allow identical calls, suffix calls, and prepending outer extensions.
- reject incompatible chains with a clear restart-session message.
- return the maximal active chain invisibly.

Implementation sketch:

```r
registerExtensions <- function(extensions, install = FALSE) {
    extensions <- validateExtensionsVector(extensions)

    old <- getRegisteredExtensions()
    relation <- compareExtensionChains(old, extensions)

    if (relation == "identical" || relation == "new_is_suffix") {
        ensureExtensionNamespaces(extensions, install = install)
        return(invisible(old))
    }

    if (relation == "old_is_suffix" || length(old) == 0) {
        ensureExtensionNamespaces(extensions, install = install)
        defineExtensionClasses(extensions)
        setRegisteredExtensions(extensions)
        return(invisible(extensions))
    }

    stop(
        "A different extension chain is already active in this session. ",
        "Please restart R before registering this chain."
    )
}
```

### `getRegisteredExtensions()`

Exported or internal helper returning the maximal active chain:

```r
getRegisteredExtensions <- function() {
    .mizerSession$extensions %||% character()
}
```

Extension packages can use this to copy the active chain into
`params@extensions`.

### `coerceToExtensionClass()`

Export this helper for extension-package constructors.

Suggested API:

```r
coerceToExtensionClass <- function(object, extensions = objectExtensions(object)) {
}
```

Behaviour:

- for `MizerParams`, coerce to the first class in `params@extensions`.
- for `MizerSim`, coerce according to `sim@params@extensions`.
- for an empty chain, coerce to the base class.
- require the object chain to be valid for the active session.

For a suffix subchain, coercion goes to the top class of that subchain, not the
top class of the maximal session chain.

### `assertExtensionChain()`

Internal helper called by public entry points after old-object upgrades.

Rules for `MizerParams`:

- `params@extensions` must be empty or a suffix of the active session chain.
- if `params@extensions` is non-empty and no compatible session chain is active,
  error and ask the user to call `registerExtensions()` or use `readParams()`.
- the object must inherit from the top class of its own `@extensions` chain.
- base objects with empty `@extensions` are valid in any session.

Rules for `MizerSim`:

- apply the same checks to `sim@params`.
- if `sim@params@extensions` is non-empty, the sim object must inherit from the
  corresponding sim-side top class.
- base `MizerSim` objects whose params have empty `@extensions` are valid in any
  session.

Do not put session-chain checks in the S4 validity functions. S4 validity
should continue to check structural object validity only. Session checks belong
in `validParams()`, `validSim()`, `readParams()`, `project()`, and other public
entry points.

## Object Lifecycle

### Creating Params Objects

Core mizer constructors can keep returning plain `MizerParams` objects.

Extension packages should follow this pattern:

```r
params <- newMultispeciesParams(...)
params@extensions <- getRegisteredExtensions()
params <- coerceToExtensionClass(params)
```

If an extension package creates params for a suffix subchain, it should set
`params@extensions` to that suffix explicitly and then call
`coerceToExtensionClass()`.

### Creating Sim Objects

`MizerSim()` and `project()` should construct the base sim object first and then
coerce it according to `sim@params@extensions`.

```r
sim <- methods::new("MizerSim", ...)
sim <- coerceToExtensionClass(sim)
```

This keeps sim dispatch aligned with the params object that produced the
simulation.

### Saving Params Objects

`saveParams()` should:

1. call `validParams()`.
2. check required extension packages and custom functions as it does now.
3. coerce the object to plain `MizerParams`.
4. preserve `params@extensions`.
5. call `saveRDS()`.

This avoids saved files depending on dynamic S4 subclass definitions.

### Reading Params Objects

`readParams()` should:

1. call `readRDS()` without immediately calling `validParams()`.
2. upgrade old objects if needed.
3. inspect `params@extensions`.
4. if the chain is non-empty, call `registerExtensions(params@extensions)`.
5. require the saved chain to be compatible with the active session chain.
6. coerce the object to the top class of its own chain.
7. call `validParams()`.
8. return the coerced object.

`readParams()` should not silently install packages by default. If installation
is desired, add an explicit argument such as:

```r
readParams(file, install_extensions = FALSE)
```

and pass it through to `registerExtensions(..., install = install_extensions)`.

### Saving And Reading Sim Objects

If mizer adds first-class sim save/read helpers later, they should mirror the
params lifecycle:

- save as plain `MizerSim`.
- preserve `sim@params@extensions`.
- read, register from `sim@params@extensions`, coerce params, coerce sim, then
  validate.

## Dispatch Refactor

### Principle

Every extension point should become an S3 generic whose first argument is the
`MizerParams` or `MizerSim` object.

Base mizer provides methods on `MizerParams` and `MizerSim`. Extension packages
provide methods on their extension classes and call `NextMethod()` to compose.

### Encounter Example

```r
projectEncounter <- function(params, n, n_pp, n_other, t = 0, ...) {
    UseMethod("projectEncounter")
}

projectEncounter.MizerParams <- function(params, n, n_pp, n_other, t = 0, ...) {
    mizerEncounter(params, n = n, n_pp = n_pp, n_other = n_other, t = t, ...)
}

projectEncounter.mizerExtA <- function(params, n, n_pp, n_other, t = 0, ...) {
    encounter <- NextMethod()
    encounter + extraEncounter(params, n, n_pp, n_other, t)
}

projectEncounter.mizerExtB <- function(params, n, n_pp, n_other, t = 0, ...) {
    encounter <- NextMethod()
    encounter * seasonalMultiplier(params, t)
}
```

With:

```r
params@extensions <- c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")
```

dispatch is:

```r
mizerExtB -> mizerExtA -> MizerParams
```

### Rates Orchestrator

Keep orchestration in core mizer:

```r
projectRates <- function(params, n, n_pp, n_other, t = 0, effort, ...) {
    UseMethod("projectRates")
}

projectRates.MizerParams <- function(params, n, n_pp, n_other,
                                     t = 0, effort, ...) {
    r <- list()
    r$encounter <- projectEncounter(params, n, n_pp, n_other, t = t, ...)
    r$feeding_level <- projectFeedingLevel(
        params, n, n_pp, n_other,
        t = t, encounter = r$encounter, ...
    )
    r$e <- projectEReproAndGrowth(
        params, n, n_pp, n_other,
        t = t,
        encounter = r$encounter,
        feeding_level = r$feeding_level,
        ...
    )
    r
}
```

Extensions may override `projectRates()`, but developer documentation should
recommend overriding lower-level generics where possible.

### First Generics To Introduce

Start with the simulation hot path:

- `projectRates`
- `projectEncounter`
- `projectFeedingLevel`
- `projectEReproAndGrowth`
- `projectERepro`
- `projectEGrowth`
- `projectPredRate`
- `projectPredMort`
- `projectFMort`
- `projectMort`
- `projectDiffusion`
- `projectRDI`
- `projectRDD`
- `projectResourceMort`
- `projectResource`
- hooks for other component dynamics where the API is stable

Then switch user-facing getters to call those same generics:

- `getRates`
- `getEncounter`
- `getFeedingLevel`
- `getPredRate`
- `getPredMort`
- `getFMort`
- `getMort`
- `getERepro`
- `getEGrowth`
- `getRDI`
- `getRDD`
- `getResourceMort`

This keeps simulation and diagnostic calculations consistent.

## Interaction With `rates_funcs`

### Phase 1

- keep `params@rates_funcs`.
- add extension-chain infrastructure.
- add S3 generics for one pilot rate, preferably encounter.
- base S3 methods may delegate to the function named in `rates_funcs`.

### Phase 2

- convert the full `projectRates()` pipeline.
- switch `project_simple()` to call `projectRates()` instead of
  `rates_fns$Rates()`.
- keep `rates_funcs` as a compatibility layer for single-function overrides.

### Phase 3

- switch user-facing getters to the S3 path.
- document `setRateFunction()` as legacy for simple, non-composable overrides.
- decide later whether `rates_funcs` remains permanently as an escape hatch.

If both mechanisms are active during transition, the S3 path is authoritative in
the simulation loop. Base S3 methods can still honour `rates_funcs`.

## Files To Change

### New Files

- `R/registerExtensions.R`
  - session-state environment.
  - `registerExtensions()`.
  - `getRegisteredExtensions()`.
  - `coerceToExtensionClass()`.
  - class-name helpers.
  - chain comparison helpers.
  - validation helpers.

- `tests/testthat/test-registerExtensions.R`
  - vector validation.
  - class creation.
  - identical registration.
  - suffix registration.
  - prepending outer extensions.
  - incompatible-chain errors.

- `tests/testthat/test-extension-dispatch.R`
  - S3 dispatch order.
  - `NextMethod()` chaining.
  - params subchain coercion.
  - sim subchain coercion via `sim@params@extensions`.

### Existing Files

- `DESCRIPTION`
  - add `R/registerExtensions.R` to `Collate:`.

- `NAMESPACE`
  - regenerate with `devtools::document()`.
  - do not edit directly.

- `R/MizerParams-class.R`
  - document `@extensions` order and suffix-subchain semantics.
  - keep S4 validity structural.
  - call extension-chain checks from `validParams.MizerParams()`.

- `R/MizerSim-class.R`
  - document that sim extension state lives in `sim@params@extensions`.
  - coerce constructed sims to the correct sim subclass.
  - call extension-chain checks from `validSim.MizerSim()`.

- `R/saveParams.R`
  - save base-class params.
  - read before validation, register from `@extensions`, coerce, then validate.

- `R/project.R`
  - ensure `project()` validates extension compatibility early.
  - ensure created `MizerSim` objects are coerced from `sim@params@extensions`.
  - later replace `rates_fns` hot-path setup with S3 dispatch.

- `R/project_methods.R`
  - convert base rate implementations into S3 methods.
  - preserve existing computational kernels where possible.

- `R/rate_functions.R`
  - rework getters to call project-rate generics.
  - preserve return classes such as `ArraySpeciesBySize`.

- `R/extension.R`
  - document the new developer-facing extension model.
  - reposition `setRateFunction()` as legacy or simple override API.

- `R/upgrade.R`
  - preserve `@extensions` during old-object upgrades.
  - avoid coercing extension objects to base except in save paths.

- `NEWS.md`
  - describe the new extension chain and migration path.

## Validation Details

### Extension Vector

`validateExtensionsVector()` should check:

- the value is a character vector.
- names are present, non-empty, unique strings.
- names are syntactically valid S4 class names.
- values are length-one strings or `NA_character_`.
- values for common extension names match exactly when comparing chains.

### Chain Comparison

Implement helpers:

```r
isSuffixChain <- function(candidate, chain) {
}

compareExtensionChains <- function(old, new) {
}
```

Comparison must use both names and values. The same extension name with a
different version/source is incompatible.

### Class Definitions

`defineExtensionClasses(extensions)` should:

- build classes from right to left.
- define params and sim classes in parallel.
- define marker subclasses only, with no additional slots.
- if a class already exists, verify it contains the expected parent.
- error immediately if a class exists with an incompatible definition.

Pseudo-code:

```r
defineExtensionClasses <- function(extensions) {
    parent_params <- "MizerParams"
    parent_sim <- "MizerSim"

    for (extension in rev(names(extensions))) {
        defineOrCheckClass(extension, parent_params)
        parent_params <- extension

        sim_class <- simExtensionClass(extension)
        defineOrCheckClass(sim_class, parent_sim)
        parent_sim <- sim_class
    }
}
```

### Object Checks

For `MizerParams`:

- `params@extensions` must be empty or a suffix of the maximal session chain.
- class must be `MizerParams` for an empty chain.
- class must inherit from `names(params@extensions)[1]` for a non-empty chain.

For `MizerSim`:

- inspect `sim@params@extensions`.
- class must be `MizerSim` for an empty params chain.
- class must inherit from `simExtensionClass(names(sim@params@extensions)[1])`
  for a non-empty params chain.

## Testing Strategy

Use unique fake extension class names in every test file, for example by adding
a random or process-specific suffix. S4 class definitions persist for the R
session and cannot be cleaned up reliably.

### Registration Tests

- empty registration leaves the chain empty and defines no classes.
- non-empty registration defines the expected params hierarchy.
- non-empty registration defines the expected sim hierarchy.
- repeated identical registration is a no-op.
- registering a suffix of the active chain is a no-op.
- registering a superchain by prepending outer extensions updates the maximal
  chain.
- registering an incompatible chain errors.
- changing the value for an already-known extension errors.

### Dispatch Tests

- fake base method returns `"base"`.
- fake extension A appends `"A"` around `NextMethod()`.
- fake extension B appends `"B"` around `NextMethod()`.
- full chain dispatch follows vector order.
- suffix object dispatch starts at the suffix top class.
- base object dispatch reaches only the base method.

### Lifecycle Tests

- extension-aware params are coerced to the top class of their own chain.
- extension-aware sims are coerced according to `sim@params@extensions`.
- `saveParams()` writes a base `MizerParams` object with preserved
  `@extensions`.
- `readParams()` registers and coerces saved extension params.
- `readParams()` accepts suffix chains in an already-registered session.
- `readParams()` rejects incompatible chains.

### Backwards Compatibility Tests

- existing models with empty `@extensions` behave unchanged.
- old saved params read correctly.
- existing `rates_funcs` overrides still work during the transition.
- tests that assign `params@rates_funcs` keep passing until the explicit
  deprecation phase.

## Implementation Order

1. Add `R/registerExtensions.R` and its tests.
2. Add `R/registerExtensions.R` to `DESCRIPTION` `Collate:`.
3. Update `MizerParams` and `MizerSim` documentation for extension semantics.
4. Implement params and sim coercion helpers.
5. Add extension-chain checks to `validParams()` and `validSim()`.
6. Update `saveParams()` and `readParams()` lifecycle.
7. Add a pilot S3 rate generic for encounter.
8. Convert the full simulation rate pipeline.
9. Switch `project_simple()` and `project()` to the new path.
10. Switch user-facing getters to the same generic path.
11. Update extension developer docs and `NEWS.md`.
12. Run `devtools::document()`, targeted tests, then broader tests.

## Immediate First Milestone

Implement infrastructure only:

- `registerExtensions()`
- chain comparison
- class definition
- params/sim coercion
- `saveParams()` and `readParams()` support
- tests proving full-chain, suffix-chain, and base-object lifecycles

Do not refactor rate functions in the first milestone. The lifecycle and class
semantics should be locked down before changing the simulation hot path.
