# Updating an extension package to the simpler S3 extension mechanism

Use this skill when updating an existing mizer extension package (such as `mizerShelf`, `mizerMR`, `mizerStarvation`, etc.) to work with mizer's S3 architecture and simplified extension dispatch mechanism.

## Overview of changes in mizer

Starting with mizer 3.3.0.9000 / 3.4:

1. **`MizerParams` and `MizerSim` are S3 classes.** They are stored as lists with class vectors `c("MizerParams")` and `c("MizerSim")`.
2. **S3 Extension Classes:** Extension packages prepend their class name to the S3 class vector (e.g. `c("mizerShelf", "MizerParams")` for params, and `c("mizerShelfSim", "MizerSim")` for simulations).
3. **Session registration is obsolete:** `mizer::registerExtension()` and dynamic S4 marker classes created at package load time have been removed.
4. **No `.onLoad` hook or active bindings needed:** Bundled datasets saved as S3 objects retain their class vectors and extension metadata under standard R lazy-loading without any active binding or `.onLoad` hook.
5. **Extension tracking via `recordExtension()`:** The extension requirement and package version are recorded in `params$extensions` (or `params@extensions`) using `recordExtension()`.
6. **Class coercion via `coerceToExtensionClass()`:** `coerceToExtensionClass(params)` reads `params$extensions` and sets the S3 class vector accordingly. `project()` automatically coerces `MizerSim` output to the corresponding `...Sim` extension class.

---

## Step-by-step update guide

### 1. Clean up `R/<pkgname>-package.R`

Remove `.onLoad()` entirely (or remove the `registerExtension()` and `makeActiveBinding()` calls within it).

- **Remove `mizer::registerExtension()`**: Session registration is no longer used.
- **Remove `makeActiveBinding()`**: Datasets in `data/` now carry their S3 class directly.
- **Remove `@importFrom methods is`**: Replace `is(x, "MizerParams")` in code with `inherits(x, "MizerParams")`.

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

### 2. Update class documentation in `R/<pkgname>-class.R`

Update the roxygen documentation for the extension class:
- Explain that the package defines **S3 extension classes** for `MizerParams` and `MizerSim`.
- Note that no `setClass()` is required.
- Mention that `new<Extension>Params()` records the extension with `recordExtension()` and calls `coerceToExtensionClass()`.

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

### 3. Update constructors and setup functions

Ensure all constructors / setup functions:
1. Call `recordExtension()` with the package name, `version`, and `requirement`.
2. Call `coerceToExtensionClass(params)` (for dispatching extensions).

```r
# In constructor:
params <- mizer::recordExtension(
    params, "myExtension",
    version = as.character(utils::packageVersion("myExtension")),
    requirement = "sizespectrum/myExtension"
)
params <- mizer::coerceToExtensionClass(params)
```

For **metadata-only** extensions (like `mizerStarvation`), call `recordExtension()` and omit `coerceToExtensionClass()`.

### 4. Migrate bundled data objects in `data/*.rda`

Bundled `MizerParams` and `MizerSim` objects in `data/` that were previously stored as S4 objects must be upgraded to S3.

Run this script once in R to upgrade and re-save the data:

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

Update `R/data.R` documentation:
```r
#' @format A [myExtension-class] object.
#'
#' The object was created with [newMyExtensionParams()] and is stored in
#' the package's `data/` directory. R's standard lazy-loading preserves its S3
#' class vector and extension metadata, so no load hook or active binding is
#' needed.
```

### 5. Replace S4 assertions and type checks

- Replace `is(x, "MizerParams")` or `is(x, "MizerSim")` with `inherits(x, "MizerParams")` or `inherits(x, "MizerSim")`.
- In `tests/testthat/`:
  - Replace `expect_s4_class(p, "MizerParams")` with `expect_s3_class(p, "MizerParams")`.
  - Replace `expect_s4_class(p, "myExtension")` with `expect_s3_class(p, "myExtension")`.
  - Replace `expect_s4_class(sim, "myExtensionSim")` with `expect_s3_class(sim, "myExtensionSim")`.

### 6. Update `DESCRIPTION`

- Remove `methods` from `Imports:` if it was only used for S4 checks / `is()`.
- Ensure `Depends:` or `Imports:` requires `mizer (>= 3.3.0.9000)` (or `>= 3.4.0`).

### 7. Update vignettes, narrative docs, and `NEWS.md`

- In vignettes (e.g. `vignettes/extension_mechanism.Rmd`), replace descriptions of S4 marker classes and `.onLoad` / `registerExtension()` with S3 extension classes and `recordExtension()` / `coerceToExtensionClass()`.
- Check vignette chunks for `plotYield(..., return_data = TRUE)` columns if `rbind` is used (the schema is `Year`, `Yield`, `Species`).
- Record the changes in `NEWS.md`:
  ```markdown
  - Replaced obsolete session extension registration and dynamic S4 marker classes with mizer's simpler S3 extension mechanism. Stored parameter objects are now S3 objects with extension class vector `c("<pkg>", "MizerParams")`, removing the need for a load hook or active binding.
  ```

### 8. Document, Test, and Check

```r
devtools::document()   # updates NAMESPACE and man/*.Rd
devtools::test()       # run test suite
devtools::check()      # verify full package build and vignettes
```
