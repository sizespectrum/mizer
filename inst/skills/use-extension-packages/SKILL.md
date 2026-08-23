---
name: use-extension-packages
description: >-
  Work with a mizer model that uses extension packages — mizerExperimental,
  mizerEcopath, mizerShelf, therMizer, mizerMR, mizerSeasonal and the rest. Use
  whenever a project loads such a package, whenever a MizerParams or MizerSim
  object has to be written to disk or read back (saveParams/readParams,
  saveSim/readSim), whenever a saved model has to be shared with or received
  from a collaborator, or whenever an extension's methods are not being called,
  a model fails to load for want of a package, or two extensions have to be
  combined and the library() order matters. Never write a params object with
  saveRDS(). To write an extension rather than use one, see the extend-mizer
  skill.
---

# Using mizer extension packages

Extension packages add new biology to mizer — extra mortality terms, additional
components such as detritus and carrion, overridden plotting functions. You do
not need to know how they are built to use them, but you do need to know how
mizer keeps track of them, because the bookkeeping is invisible until it goes
wrong.

Two rules cover almost everything:

1. **Load the extension packages with `library()` before you build, read or use
   a model that needs them.** The chain is rebuilt in load order every session.
2. **Persist models with `saveParams()`/`readParams()` (or
   `saveSim()`/`readSim()`), never with `saveRDS()`/`readRDS()`.** The file on
   disk deliberately holds a *base-class* object; the extension class is
   restored on the way back in by `readParams()`. So a bare `readRDS()` of a
   perfectly good file returns an object whose extension methods are no longer
   called — the model quietly gives base-mizer answers instead of failing.

## Available extension packages

- **[mizerExperimental](https://github.com/sizespectrum/mizerExperimental)** —
  A community-driven collection of experimental features being refined before
  potential inclusion in the core mizer package.
- **[mizerShiny](https://github.com/gustavdelius/mizerShiny)** — Provides a
  web-based Shiny interface for running and exploring mizer models without
  writing R code.
- **[therMizer](https://github.com/sizespectrum/therMizer)** — Incorporates
  temperature-dependent metabolic rates and aerobic scope into mizer, with
  support for size- and depth-varying thermal exposure.
- **[mizerEcopath](https://github.com/gustavdelius/mizerEcopath)** — Provides
  tools for calibrating mizer models using Ecopath biomass estimates, bridging
  the two modelling frameworks.
- **[mizerStarvation](https://github.com/sizespectrum/mizerStarvation)** —
  Implements starvation mortality, allowing fish to die when food availability
  is insufficient to sustain their metabolic needs.
- **[mizerMR](https://github.com/sizespectrum/mizerMR)** — Supports multiple
  size-structured background resource spectra, enabling species to have
  different prey preferences and ontogenetic dietary shifts.
- **[mizerShelf](https://github.com/sizespectrum/mizerShelf)** — Adds detritus
  and carrion components for more realistic benthic ecosystem modelling on
  continental shelves.
- **[mizerEvolution](https://github.com/baldrech/mizerEvolution)** — Enables
  simulation of evolutionary trait changes and species invasions by treating
  species as pools of phenotypes subject to natural selection.
- **[mizerSeasonal](https://github.com/gustavdelius/mizerSeasonal)** —
  Introduces seasonal dynamics by simulating gonadal mass accumulation and
  spawning events over the course of a year.
- **[mizerStomach](https://github.com/gustavdelius/mizerStomach)** — Helps
  determine predator size-selectivity functions by fitting them to stomach
  content data.

## The extension chain

Every mizer extension package announces itself to mizer when it is loaded.
Internally, mizer maintains a list of all the currently loaded extensions, in
the order they registered themselves. This list is called the **extension
chain**.

```r
getRegisteredExtensions()
```

The result is a named character vector. The **names** are the extension
identifiers (which double as S4 class names). The **values** are installation
specifications — a version string such as `"1.2.0"` for packages on CRAN, or a
GitHub path such as `"sizespectrum/mizerShelf"` for packages that are only on
GitHub. The first element is the *outermost* extension (the one with the highest
dispatch priority), the last is the innermost.

If no extension packages are loaded, `getRegisteredExtensions()` returns an
empty, unnamed character vector.

## Why loading order matters

When two extensions both override the same mizer function — say `getBiomass()` —
R decides which version to call first from the class of the `params` object,
working through the class hierarchy from outermost to innermost. The outermost
extension gets the first say, and each extension can modify the result of the
extensions below it.

So **the package loaded last is outermost**:

```r
library(mizerShelf)   # innermost
library(mizerFoo)     # outermost: its getBiomass() runs first and wraps the rest
```

Both extensions contribute either way, so in most cases the two orderings give
the same answer, because the extensions touch different parts of the
calculation. The order matters only when the extensions interact — when one
extension's contribution depends on what another computed. Check each package's
documentation for a required load order.

## Working with the extension chain

### Clearing the registry

To replace the chain explicitly without restarting R, first clear the session's
registry:

```r
clearExtensionChain()
```

`getRegisteredExtensions()` then returns an empty vector. Clearing the chain does
**not** unload the packages themselves; it only removes their entries from
mizer's registry. Calling `library()` for a package that is still loaded does
not reliably rerun its `.onLoad()` hook, so either restart or unload and load the
packages again, or restore the chain explicitly with `registerExtensions()`
before creating or using dependent params objects.

### Setting the chain manually

To take full control — for example to reproduce the exact chain recorded in a
collaborator's params object — set it directly with `registerExtensions()`:

```r
chain <- c(mizerFoo = "owner/mizerFoo", mizerShelf = "sizespectrum/mizerShelf")
registerExtensions(chain)
```

The names must match the identifiers each package uses when it registers itself
(normally the package name). The values are installation specifications, used
when a missing package has to be installed automatically. This is an advanced
manoeuvre: in normal use you call `library()` for each extension and the chain
builds itself.

## Params objects and the extension chain

When an extension package creates a `MizerParams` object (for example
`mizerShelf::newDetritusCarrionParams()`), it records the extension packages
actually applied to that object in the `extensions` slot:

```r
params@extensions
```

That record serves two purposes:

1. **Reproducibility.** It says which extension packages built the model, the
   installation requirement for each, and the package version whose object
   layout each component conforms to.
2. **Class restoration.** On reload, mizer uses it to promote the object back to
   the correct S4 class, so that functions like `getBiomass()` keep dispatching
   to the right extension methods.

### Saving and loading

```r
saveParams(params, "my_model.rds")
params <- readParams("my_model.rds")
```

Before saving, `saveParams()` warns if the model relies on custom functions
defined only in your R session — a custom rate function, selectivity function or
predation kernel written in a script. Those functions are not stored in the
file, so the script has to travel with the `.rds`.

`readParams()` upgrades the object if it was written by an older mizer, reads
`params@extensions`, calls `registerExtensions()` so the session knows the saved
chain, and promotes the object to the correct S4 class. As long as the required
packages are installed, this is seamless; if one is missing, `readParams()` stops
with an error naming it.

`MizerSim` objects carry their params object inside them, so the chain is
embedded there too. Use `saveSim()` and `readSim()`, which do the same
registration and coercion:

```r
sim <- project(params, t_max = 10)
saveSim(sim, "my_simulation.rds")
sim <- readSim("my_simulation.rds")
```

For the metadata that should accompany a model you intend to share
(`setMetadata()`), see the `build-model` skill.

## Sharing models with collaborators

A collaborator needs the same extension packages installed. `readParams()` and
`readSim()` tell them which are missing. To install them automatically from the
specifications stored in `params@extensions`:

```r
params <- readParams("my_model.rds", install_extensions = TRUE)
```

This installs each package from the recorded requirement specification, such as
a CRAN version requirement or a GitHub repository. The separate recorded
version stamp describes the object layout for upgrade purposes; it does not
necessarily pin installation to that exact package release.

## Built-in example models from extension packages

Extension packages often ship ready-made `MizerParams` or `MizerSim` objects as
example models. As long as the package follows mizer's conventions, these work
as soon as the package is loaded:

```r
library(mizerShelf)
NWMed_params   # already has the correct extension class -- no extra steps needed
```

An object from an older package that does not behave as expected can be repaired
with `coerceToExtensionClass()`.

## Common scenarios

### A model fails to load for want of a package

```
Error in readParams("my_model.rds") :
  Some required extension packages are not installed: mizerShelf
```

Install it manually (`pak::pkg_install("sizespectrum/mizerShelf")`) or reload
with `readParams("my_model.rds", install_extensions = TRUE)`.

### Checking what a params object requires

```r
params@extensions
```

The names are the extension identifiers. Current entries each contain a
`requirement` (where or at what minimum version to install the package) and a
`version` stamp (the package version whose object layout the component conforms
to). Legacy objects may still show the older named-character-vector form;
mizer accepts both.

### Packages were loaded in the wrong order

Either restart R and load the packages in the desired order, or explicitly
register the desired chain. Merely calling `library()` again after
`clearExtensionChain()` is insufficient when those namespaces are still
loaded. `registerExtensions()` takes the chain in outermost-first order:

```r
clearExtensionChain()
registerExtensions(c(
    mizerFoo = "owner/mizerFoo",                  # outermost
    mizerShelf = "sizespectrum/mizerShelf"        # innermost
))
```

Afterwards, create or re-read the params objects that use that chain. In a fresh
session the equivalent loading order is `library(mizerShelf)` followed by
`library(mizerFoo)`.

### Making a script reproducible

Put the `library()` calls at the top in a fixed order. Each package registers
itself in `.onLoad`, so sourcing the script fresh always rebuilds the same chain.

<!-- agent-only -->
## Symptom index

| Symptom | Cause | Fix |
|---|---|---|
| `readParams()`/`readSim()` errors "Some required extension packages are not installed" | The chain recorded in the file names a package this library does not have | Install it, or re-read with `install_extensions = TRUE` |
| An extension's method is not called: results match plain mizer, no error | The object lost its marker class — usually read with `readRDS()` instead of `readParams()`, or the package was not loaded | `library(<pkg>)`, then `params <- coerceToExtensionClass(params)`; re-read the file with `readParams()` in future |
| `params@extensions` is empty on a model an extension built | The extension's setup function did not record itself, or the object predates that mechanism | Load the package and rerun its setup/conversion function or rebuild the model; coercion alone does not persist the missing record |
| `saveParams()` warns "Your model is using the functions …" | A custom rate, selectivity or kernel function lives only in the session | Ship the defining script alongside the `.rds`; the warning is not a failure |
| Results change depending on the order of the `library()` calls | Two extensions override the same function and interact | Fix the order deliberately in a fresh session, or call `registerExtensions()` with an explicit outermost-first chain, then rebuild or re-read the objects |
| `getRegisteredExtensions()` is empty although the package is attached | The registry was cleared after the namespace loaded | Restart or unload and load the package so `.onLoad()` runs, or call `registerExtensions()` explicitly |

Before writing a params object to disk in a project that loads any extension
package, check that the code path uses `saveParams()`. `saveRDS()` produces a
file that loads without error and silently dispatches to base mizer.
<!-- /agent-only -->

## See also

- The `extend-mizer` skill — the mechanisms an extension is built from.
- [Creating a mizer extension package](guide-create-extension-package.html) — how to
  package your own extension and make it composable with others.
- `?registerExtension`, `?registerExtensions`, `?getRegisteredExtensions`,
  `?clearExtensionChain`, `?coerceToExtensionClass`
- `?saveParams`, `?readParams`, `?saveSim`, `?readSim`
