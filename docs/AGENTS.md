# AI Agent Instructions for mizer

mizer is an R package for dynamic multi-species size-spectrum modelling
of fish communities.

## Architecture

**`MizerParams`** (S4) — central object passed to nearly all functions.
Modified via setter functions that return new copies:
`setFishing(params, ...)`, etc.

**`MizerSim`** (S4) — stores simulation output: arrays `n` (time ×
species × size), `n_pp`, `n_other`, `effort`, and the `MizerParams`
used.

**`ArraySpeciesBySize`** / **`ArrayTimeBySpecies`** /
**`ArrayTimeBySpeciesBySize`** (S3) — wrap some arrays with metadata.
When assigning back to S4 slots, use `slot[] <- value` (not
`slot <- value`) to strip the class.

**`species_params`** / **`given_species_params`** / **`gear_params`**
(S3) — parameter tables subclassing `data.frame` stored in S4 slots.
Users should access/modify these tables via S3 generics
(e.g. `species_params(params)` or `gear_params(params)`), but package
code and developers can modify slots directly
(e.g. `params@species_params`) when appropriate. Inline modifications on
the S3 objects (via `[<-`, `$<-`, `[[<-` S3 methods) preserve the
subclass and trigger reactive validation checks and conversions
(e.g. length-to-weight). `given_species_params` holds what the user
supplied **plus the defaults of any function argument that sets a
species parameter** (e.g. `n` and `p` from
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)),
even when the user did not override the argument; defaults that are not
function arguments stay out of it.

**Customisable rate functions**: users replace rate functions by storing
a custom function name in `params@rates_funcs`. Dispatch via
`get(params@rates_funcs$FunctionName)(params, ...)`.

**Auto-generated files** — never edit `NAMESPACE`, `man/`,
`RcppExports.R`, or `RcppExports.cpp` directly.

## Code Conventions

- **Indentation**: 4 spaces
- **Naming**: camelCase or snake_case for functions/variables;
  PascalCase for classes
- **Language**: British English (en-GB) — “colour”, “behaviour”,
  “modelling”
- When documenting a mizer S3 generic whose methods share a man page
  (combined with `@rdname`/`@name`), follow the steps in
  `.claude/skills/document-s3-generics.md`.
- When adding, moving or removing a species parameter default, follow
  `.claude/skills/species-param-defaults.md`. A default belongs to the
  rate setter that reads the parameter; only parameters that no single
  rate setter owns are defaulted centrally.

## Testing

- Use `expect_doppelganger()` (vdiffr) for plot tests

- Use snapshot tests for complex outputs

- Run
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  after adding or changing exports

- Run
  [`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
  before running tests

- Run only relevant tests with `devtools::test(filter = "pattern")`.
  Running all tests is too slow.

- Lint a file with
  [`lintr::lint()`](https://lintr.r-lib.org/reference/lint.html)

- After editing C++ source:
  [`devtools::clean_dll(); devtools::load_all()`](https://pkgbuild.r-lib.org/reference/clean_dll.html)

- After modifying the `MizerParams` or `MizerSim` class (new/removed
  slots, changes to `@rates_funcs`, etc.), follow the steps in
  `.claude/skills/upgrade-mizer-data.md`.

- When snapshot tests fail because values legitimately changed, run
  [`testthat::snapshot_accept()`](https://testthat.r-lib.org/reference/snapshot_accept.html)
  to promote the `.new.md` files into the canonical `.md` snapshots.

- Avoid top-level [`data()`](https://rdrr.io/r/utils/data.html) calls in
  test files. `testthat 3.x` shares `.GlobalEnv` across all tests. Use
  the `_small` fixtures from `helper.R` instead.

### Test helper objects

`helper.R` defines small objects used across all test files. All are
named with the `_small` suffix to avoid shadowing the built-in package
datasets:

- **`NS_params_small`**: 3 species (Sprat=1, Herring=2, Cod=3), 20 size
  bins, 3 gears (Industrial effort=0, Pelagic effort=1, Otter
  effort=0.5). Sprat is always unfished (Industrial gear, zero effort).
- **`NS_sim_small`**: `project(NS_params_small, t_max=3, t_save=1)` — 4
  time steps at t=0,1,2,3.
- **`NS_species_params_small`**, **`NS_species_params_gears_small`**,
  **`inter_small`**: 3-species subsets of the full package datasets
  (rows 1, 4, 11 = Sprat, Herring, Cod).

Tests that hardcode species indices (e.g. `[12, ]`) or time indices
(e.g. `times[10]`) must use values valid for 3 species and 4 time steps.

## Before Submitting

- After adding a new file under `R/`, add it to the `Collate:` field in
  `DESCRIPTION` (roxygen2 does not manage this automatically in this
  package).
- After adding or renaming exported functions, add them to the
  `appropriate` section in `pkgdown/_pkgdown.yml` so they appear on the
  reference page of the website.
- Update `NEWS.md` when adding features or fixing bugs.
- When creating a pull request, always include the summary of your
  changes in the PR body.
