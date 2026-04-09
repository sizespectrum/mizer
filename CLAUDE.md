# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About mizer

mizer is an R package for dynamic multi-species size-spectrum modelling of fish communities. It models marine ecosystems subject to fishing, tracking individual fish growth from egg size to maximum size and capturing ontogenetic diet shifts.

## Common Commands

```r
devtools::load_all()        # Load package for development
devtools::document()        # Regenerate NAMESPACE and man/ pages from roxygen2
devtools::test()            # Run all tests
devtools::check()           # Full R CMD check
lintr::lint_package()       # Lint the package

# Run a single test file
testthat::test_file("tests/testthat/test-filename.R")

# After editing C++ source
devtools::clean_dll(); devtools::load_all()
```

## Architecture

### Core Classes

**`MizerParams`** (S4, `R/MizerParams-class.R`) — the central object passed to nearly all functions. Holds all model configuration: species parameters, size grids (`w`, `w_full`), interaction matrices, gear selectivity, rate function overrides (`@rates_funcs`), and resource dynamics. Validated by `validMizerParams()`. Modified via setter functions that return new copies: `setFishing(params, ...)`, `setInteraction(params, ...)`, etc.

**`MizerSim`** (S4, `R/MizerSim-class.R`) — stores simulation output: a 3D array `n` (time × species × size), `n_pp` (time × size) for resource, `n_other` for additional biomass components, `effort` history, and the `MizerParams` used.

**`ArraySpeciesBySize`** (S3, `R/ArraySpeciesBySize-class.R`) — wraps 2D arrays (species × size) returned by mizer functions with metadata (`value_name`, `units`). Inherits from `matrix`/`array`. Provides enhanced `print()`, `summary()`, `plot()`, and `as.data.frame()`.

### Execution Flow

1. **Setup**: `newMultispeciesParams()` / `newSingleSpeciesParams()` / `newCommunityParams()`
2. **Configure**: `set*()` functions (`setFishing()`, `setPredKernel()`, `setMetabolicRate()`, …)
3. **Rates**: `getEncounter()`, `getFeedingLevel()`, `getPredMort()`, `getFMort()`, `getRates()` — each returns an `ArraySpeciesBySize` matrix (species × size)
4. **Project**: `project(params, t_max = 100, effort = ...)` → `MizerSim`
5. **Analyse**: `getYield()`, `getBiomass()`, `getSSB()`, `plotSpectra()`, etc.

### Customisable Rate Functions

Users can replace any rate function by storing a custom function name in `params@rates_funcs`. Calls dispatch via `get(params@rates_funcs$FunctionName)(params, ...)`. This is the primary extensibility mechanism.

### C++ Integration

Performance-critical inner projection loop lives in `src/inner_project_loop.cpp` and `src/project_n_loop.cpp`. `RcppExports.R` and `RcppExports.cpp` are auto-generated — never edit them directly.

### Extensibility via "Other" Components

Arbitrary biomass components beyond species can be added (e.g., detritus, zooplankton). Stored in `n_other`; custom rate functions can be registered for them.

## Code Conventions

- **Indentation**: 4 spaces
- **Naming**: camelCase or snake_case for functions/variables; PascalCase for classes
- **Language**: British English (en-GB) — "colour", "behaviour", "modelling"
- **Documentation**: All exported functions require roxygen2 with `@param`, `@return`, `@export`; use `@seealso` for cross-references
- **Tests**: testthat edition 3; use `expect_doppelganger()` (vdiffr) for plot tests, snapshot tests for complex outputs
- After adding exports, run `devtools::document()` to regenerate `NAMESPACE`
- Update `NEWS.md` when adding features or fixing bugs
