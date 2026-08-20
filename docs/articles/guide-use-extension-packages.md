# Guide: Using mizer extension packages

This guide explains what happens when you load one or more mizer
extension packages, why the order in which you load them can matter, and
how to save and share models that use extensions. To write an extension
rather than use one, see the [guide to creating a mizer extension
package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.md).

Extension packages add new biology to mizer — extra mortality terms,
additional components such as detritus and carrion, overridden plotting
functions. You do not need to know how they are built to use them, but
you do need to know how mizer keeps track of them, because the
bookkeeping is invisible until it goes wrong.

Two rules cover almost everything:

1.  **Load the extension packages with
    [`library()`](https://rdrr.io/r/base/library.html) before you build,
    read or use a model that needs them.** The chain is rebuilt in load
    order every session.
2.  **Persist models with
    [`saveParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)/[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
    (or
    [`saveSim()`](https://sizespectrum.org/mizer/reference/saveParams.md)/[`readSim()`](https://sizespectrum.org/mizer/reference/saveParams.md)),
    never with
    [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/[`readRDS()`](https://rdrr.io/r/base/readRDS.html).**
    The file on disk deliberately holds a *base-class* object; the
    extension class is restored on the way back in by
    [`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md).
    So a bare [`readRDS()`](https://rdrr.io/r/base/readRDS.html) of a
    perfectly good file returns an object whose extension methods are no
    longer called — the model quietly gives base-mizer answers instead
    of failing.

------------------------------------------------------------------------

## Available extension packages

- **[mizerExperimental](https://github.com/sizespectrum/mizerExperimental)**
  — A community-driven collection of experimental features being refined
  before potential inclusion in the core mizer package.
- **[mizerShiny](https://github.com/gustavdelius/mizerShiny)** —
  Provides a web-based Shiny interface for running and exploring mizer
  models without writing R code.
- **[therMizer](https://github.com/sizespectrum/therMizer)** —
  Incorporates temperature-dependent metabolic rates and aerobic scope
  into mizer, with support for size- and depth-varying thermal exposure.
- **[mizerEcopath](https://github.com/gustavdelius/mizerEcopath)** —
  Provides tools for calibrating mizer models using Ecopath biomass
  estimates, bridging the two modelling frameworks.
- **[mizerStarvation](https://github.com/sizespectrum/mizerStarvation)**
  — Implements starvation mortality, allowing fish to die when food
  availability is insufficient to sustain their metabolic needs.
- **[mizerMR](https://github.com/sizespectrum/mizerMR)** — Supports
  multiple size-structured background resource spectra, enabling species
  to have different prey preferences and ontogenetic dietary shifts.
- **[mizerShelf](https://github.com/sizespectrum/mizerShelf)** — Adds
  detritus and carrion components for more realistic benthic ecosystem
  modelling on continental shelves.
- **[MizerEvolution](https://github.com/baldrech/MizerEvolution)** —
  Enables simulation of evolutionary trait changes and species invasions
  by treating species as pools of phenotypes subject to natural
  selection.
- **[mizerSeasonal](https://github.com/gustavdelius/mizerSeasonal)** —
  Introduces seasonal dynamics by simulating gonadal mass accumulation
  and spawning events over the course of a year.
- **[mizerStomach](https://github.com/gustavdelius/mizerStomach)** —
  Helps determine predator size-selectivity functions by fitting them to
  stomach content data.

------------------------------------------------------------------------

## The extension chain

Every mizer extension package announces itself to mizer when it is
loaded. Internally, mizer maintains a list of all the currently loaded
extensions, in the order they registered themselves. This list is called
the **extension chain**.

``` r

getRegisteredExtensions()
```

The result is a named character vector. The **names** are the extension
identifiers (which double as S4 class names). The **values** are
installation specifications — a version string such as `"1.2.0"` for
packages on CRAN, or a GitHub path such as `"sizespectrum/mizerShelf"`
for packages that are only on GitHub. The first element is the
*outermost* extension (the one with the highest dispatch priority), the
last is the innermost.

If no extension packages are loaded,
[`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md)
returns an empty, unnamed character vector.

------------------------------------------------------------------------

## Why loading order matters

When two extensions both override the same mizer function — say
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
— R decides which version to call first from the class of the `params`
object, working through the class hierarchy from outermost to innermost.
The outermost extension gets the first say, and each extension can
modify the result of the extensions below it.

So **the package loaded last is outermost**:

``` r

library(mizerShelf)   # innermost
library(mizerFoo)     # outermost: its getBiomass() runs first and wraps the rest
```

Both extensions contribute either way, so in most cases the two
orderings give the same answer, because the extensions touch different
parts of the calculation. The order matters only when the extensions
interact — when one extension’s contribution depends on what another
computed. Check each package’s documentation for a required load order.

------------------------------------------------------------------------

## Working with the extension chain

### Restarting with a clean chain

If the packages were loaded in the wrong order, clear the session’s
registry rather than restarting R:

``` r

clearExtensionChain()
```

[`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md)
then returns an empty vector. Clearing the chain does **not** unload the
packages themselves; it only removes their entries from mizer’s
registry, so you must reload or re-register the extensions before
creating or using params objects that depend on them.

### Setting the chain manually

To take full control — for example to reproduce the exact chain recorded
in a collaborator’s params object — set it directly with
[`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md):

``` r

chain <- c(mizerFoo = "owner/mizerFoo", mizerShelf = "sizespectrum/mizerShelf")
registerExtensions(chain)
```

The names must match the identifiers each package uses when it registers
itself (normally the package name). The values are installation
specifications, used when a missing package has to be installed
automatically. This is an advanced manoeuvre: in normal use you call
[`library()`](https://rdrr.io/r/base/library.html) for each extension
and the chain builds itself.

------------------------------------------------------------------------

## Params objects and the extension chain

When an extension package creates a
[`MizerParams`](https://sizespectrum.org/mizer/reference/MizerParams.md)
object (for example
[`mizerShelf::newDetritusCarrionParams()`](https://sizespectrum.org/mizerShelf/reference/newDetritusCarrionParams.html)),
it stamps the object with the full extension chain active at that
moment, stored in the `extensions` slot:

``` r

params@extensions
```

That record serves two purposes:

1.  **Reproducibility.** It says which extension packages, and which
    versions, built the model.
2.  **Class restoration.** On reload, mizer uses it to promote the
    object back to the correct S4 class, so that functions like
    [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
    keep dispatching to the right extension methods.

### Saving and loading

``` r

saveParams(params, "my_model.rds")
params <- readParams("my_model.rds")
```

Before saving,
[`saveParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
warns if the model relies on custom functions defined only in your R
session — a custom rate function, selectivity function or predation
kernel written in a script. Those functions are not stored in the file,
so the script has to travel with the `.rds`.

[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
upgrades the object if it was written by an older mizer, reads
`params@extensions`, calls
[`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md)
so the session knows the saved chain, and promotes the object to the
correct S4 class. As long as the required packages are installed, this
is seamless; if one is missing,
[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
stops with an error naming it.

[`MizerSim`](https://sizespectrum.org/mizer/reference/MizerSim.md)
objects carry their params object inside them, so the chain is embedded
there too. Use
[`saveSim()`](https://sizespectrum.org/mizer/reference/saveParams.md)
and
[`readSim()`](https://sizespectrum.org/mizer/reference/saveParams.md),
which do the same registration and coercion:

``` r

sim <- project(params, t_max = 10)
saveSim(sim, "my_simulation.rds")
sim <- readSim("my_simulation.rds")
```

For the metadata that should accompany a model you intend to share
([`setMetadata()`](https://sizespectrum.org/mizer/reference/setMetadata.md)),
see the [guide to building a mizer
model](https://sizespectrum.org/mizer/articles/guide-build-model.md).

------------------------------------------------------------------------

## Sharing models with collaborators

A collaborator needs the same extension packages installed.
[`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
and
[`readSim()`](https://sizespectrum.org/mizer/reference/saveParams.md)
tell them which are missing. To install them automatically from the
specifications stored in `params@extensions`:

``` r

params <- readParams("my_model.rds", install_extensions = TRUE)
```

This fetches the recorded version of each package from CRAN or GitHub as
appropriate.

------------------------------------------------------------------------

## Built-in example models from extension packages

Extension packages often ship ready-made `MizerParams` or `MizerSim`
objects as example models. As long as the package follows mizer’s
conventions, these work as soon as the package is loaded:

``` r

library(mizerShelf)
NWMed_params   # already has the correct extension class -- no extra steps needed
```

An object from an older package that does not behave as expected can be
repaired with
[`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md).

------------------------------------------------------------------------

## Common scenarios

### A model fails to load for want of a package

    Error in readParams("my_model.rds") :
      Some required extension packages are not installed: mizerShelf

Install it manually (`pak::pkg_install("sizespectrum/mizerShelf")`) or
reload with `readParams("my_model.rds", install_extensions = TRUE)`.

### Checking what a params object requires

``` r

params@extensions
```

The names are the package names; the values say where to get each one.

### Packages were loaded in the wrong order

Call
[`clearExtensionChain()`](https://sizespectrum.org/mizer/reference/clearExtensionChain.md),
then reload in the order you want — no need to restart R. Afterwards,
create or re-read your params objects so they pick up the new chain.

``` r

clearExtensionChain()
library(mizerShelf)   # innermost (loaded first)
library(mizerFoo)     # outermost (loaded last)
```

### Making a script reproducible

Put the [`library()`](https://rdrr.io/r/base/library.html) calls at the
top in a fixed order. Each package registers itself in `.onLoad`, so
sourcing the script fresh always rebuilds the same chain.

------------------------------------------------------------------------

## See also

- The [guide to extending
  mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.md)
  — the mechanisms an extension is built from.
- [Creating a mizer extension
  package](https://sizespectrum.org/mizer/articles/guide-create-extension-package.md)
  — how to package your own extension and make it composable with
  others.
- [`?registerExtension`](https://sizespectrum.org/mizer/reference/registerExtension.md),
  [`?registerExtensions`](https://sizespectrum.org/mizer/reference/registerExtensions.md),
  [`?getRegisteredExtensions`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md),
  [`?clearExtensionChain`](https://sizespectrum.org/mizer/reference/clearExtensionChain.md),
  [`?coerceToExtensionClass`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md)
- [`?saveParams`](https://sizespectrum.org/mizer/reference/saveParams.md),
  [`?readParams`](https://sizespectrum.org/mizer/reference/saveParams.md),
  [`?saveSim`](https://sizespectrum.org/mizer/reference/saveParams.md),
  [`?readSim`](https://sizespectrum.org/mizer/reference/saveParams.md)
