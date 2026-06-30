# AI Agent Instructions for mizer

mizer is an R package for dynamic multi-species size-spectrum modelling of fish communities.

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

**`MizerParams`** (S4) — central object passed to nearly all functions. Modified via setter functions that return new copies: `setFishing(params, ...)`, etc.

**`MizerSim`** (S4) — stores simulation output: arrays `n` (time × species × size), `n_pp`, `n_other`, `effort`, and the `MizerParams` used.

**`ArraySpeciesBySize`** / **`ArrayTimeBySpecies`** / **`ArrayTimeBySpeciesBySize`** (S3) — wrap some arrays with metadata. When assigning back to S4 slots, use `slot[] <- value` (not `slot <- value`) to strip the class.

**Customisable rate functions**: users replace rate functions by storing a custom function name in `params@rates_funcs`. Dispatch via `get(params@rates_funcs$FunctionName)(params, ...)`.

**Auto-generated files** — never edit `NAMESPACE`, `man/`, `RcppExports.R`, or `RcppExports.cpp` directly.

## Code Conventions

- **Indentation**: 4 spaces
- **Naming**: camelCase or snake_case for functions/variables; PascalCase for classes
- **Language**: British English (en-GB) — "colour", "behaviour", "modelling"
- When documenting a mizer S3 generic whose methods share a man page (combined
  with `@rdname`/`@name`), follow the steps in
  `.claude/skills/document-s3-generics.md`. Claude Code users can invoke this as
  `/document-s3-generics`.

## Testing

- Use `expect_doppelganger()` (vdiffr) for plot tests
- Use snapshot tests for complex outputs
- Run `devtools::document()` after adding or changing exports
- Run `devtools::load_all()` before running tests
- After modifying the `MizerParams` or `MizerSim` class (new/removed slots, changes to `@rates_funcs`, etc.), follow the steps in `.claude/skills/upgrade-mizer-data.md`. Claude Code users can invoke this as `/upgrade-mizer-data`.
- When snapshot tests fail because values legitimately changed, run `testthat::snapshot_accept()` to promote the `.new.md` files into the canonical `.md` snapshots.

### Test helper objects

`helper.R` defines small objects used across all test files. All are named with the `_small` suffix to avoid shadowing the built-in package datasets:

- **`NS_params_small`**: 3 species (Sprat=1, Herring=2, Cod=3), 20 size bins, 3 gears (Industrial effort=0, Pelagic effort=1, Otter effort=0.5). Sprat is always unfished (Industrial gear, zero effort).
- **`NS_sim_small`**: `project(NS_params_small, t_max=3, t_save=1)` — 4 time steps at t=0,1,2,3.
- **`NS_species_params_small`**, **`NS_species_params_gears_small`**, **`inter_small`**: 3-species subsets of the full package datasets (rows 1, 4, 11 = Sprat, Herring, Cod).

Tests that hardcode species indices (e.g. `[12, ]`) or time indices (e.g. `times[10]`) must use values valid for 3 species and 4 time steps.

### Global environment and `data()` in test files

In testthat 3.x all test files share the same global environment. A top-level `data()` call in one file writes to `.GlobalEnv` and therefore affects all alphabetically-later files. `test-backwards_compatibility.R` deliberately loads the full 12-species `NS_species_params_gears` and `inter` datasets at the top for its own snapshot tests — those names are intentionally **not** renamed to `_small` in that file. No restore block is needed at the bottom because the `_small` test fixtures have distinct names and are not overwritten.

## Before Submitting

- After adding a new file under `R/`, add it to the `Collate:` field in `DESCRIPTION` (roxygen2 does not manage this automatically in this package).
- After adding or renaming exported functions, add them to the `appropriate` section in `pkgdown/_pkgdown.yml` so they appear on the reference page of the website.
- Update `NEWS.md` when adding features or fixing bugs.
- When creating a pull request, always include the summary of your changes in the PR body.

