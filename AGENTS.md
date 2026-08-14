# AI Agent Instructions for mizer

mizer is an R package for dynamic multi-species size-spectrum modelling of fish communities.

## Architecture

**`MizerParams`** (S4) — central object passed to nearly all functions. Modified via setter functions that return new copies: `setFishing(params, ...)`, etc.

**`MizerSim`** (S4) — stores simulation output: arrays `n` (time × species × size), `n_pp`, `n_other`, `effort`, and the `MizerParams` used.

**`ArraySpeciesBySize`** / **`ArrayTimeBySpecies`** / **`ArrayTimeBySpeciesBySize`** (S3) — wrap some arrays with metadata. When assigning back to S4 slots, use `slot[] <- value` (not `slot <- value`) to strip the class.

**`species_params`** / **`given_species_params`** / **`gear_params`** (S3) — parameter tables subclassing `data.frame` stored in S4 slots. Users should access/modify these tables via S3 generics (e.g. `species_params(params)` or `gear_params(params)`), but package code and developers can modify slots directly (e.g. `params@species_params`) when appropriate. Editing one of these tables on its own does nothing but edit it: the subclass is preserved (by the base `data.frame` methods — there are deliberately no `[<-`, `$<-` or `[[<-` methods), and no validation, conversion or warning happens. The checks and conversions (misspelled column names, length-to-weight, consistency) run when the table is validated, which is what `species_params()`/`gear_params()` do to a plain data frame and what `species_params<-()`/`gear_params<-()` do to the table they are given. This is also what makes the length/weight precedence rule work: `species_params<-()` compares the incoming table against the model's, so the value you changed wins; a table edited on its own carries no such history, so a length and a weight that both differ count as given at the same time and the weight wins. `given_species_params` holds what the user supplied **plus the defaults of any function argument that sets a species parameter** (e.g. `n` and `p` from `newMultispeciesParams()`), even when the user did not override the argument; defaults that are not function arguments stay out of it.

**Customisable rate functions**: users replace rate functions by storing a custom function name in `params@rates_funcs`. Dispatch via `get(params@rates_funcs$FunctionName)(params, ...)`.

**Two quadrature schemes**: the `second_order_w` slot selects how the model is discretised on the size grid — `flux` picks the advective reconstruction, and `bin_average` decides whether size-dependent factors are integrated over their bin or point-sampled at the left bin boundary. Both default to the first-order scheme, so previous mizer versions reproduce byte-for-byte. This is invisible in the code you are likely to be reading: a function that ignores `bin_average` still looks correct and still passes its tests on the default path. Anything that integrates over the size grid must handle both schemes and be tested under both.

**Gated validation**: `validParams()` splits into three jobs — upgrade (gated on the version string), repair plus `validObject()` (gated on a fingerprint), and the unconditional `is.finite()` scans of the rate arrays. The fingerprint, `validation_key()`, is a hash of every slot that `repair_params()` reads or writes and every slot that `validMizerParams()` inspects, taking only the dimensions and dimension names of the big rate arrays, never their values. It is recalculated on every call and never stored on the object, so it cannot go stale. This gate stays sound only as long as that correspondence holds: **if you make `repair_params()` or `validMizerParams()` depend on a slot, add that slot to `validation_key()`**, or the new logic will be silently skipped for objects mizer has already seen. Values inside the rate arrays are deliberately out of the fingerprint and are guarded by `check_finite()` instead.


**Auto-generated files** — never edit `NAMESPACE`, `man/`, `RcppExports.R`, or `RcppExports.cpp` directly. The `vignettes/cheatsheet-*.Rmd` articles are also generated: their single source is `inst/skills/<topic>/SKILL.md`, which doubles as the agent skill. Edit the skill and re-run `source("dev_scripts/build_cheatsheets.R"); build_cheatsheets()`. Content that belongs in the skill but not in the article (an agent's diagnostic procedure, a lookup table keyed by symptom rather than by topic) goes in an `<!-- agent-only -->` block, which the generator drops. Likewise `docs/llms.txt` and `inst/llms.txt` are generated: pkgdown writes a raw `docs/llms.txt` during `build_site()`, and `source("dev_scripts/build_llms.R"); build_llms()` then swaps its README-derived preamble for `pkgdown/llms-header.md` — the file to edit — and copies the result to `inst/llms.txt`. That installed copy is how `mizerAgents::setup_mizer_agent()` points an agent at the API index for the mizer it has, so it must ship: do not add it to `.Rbuildignore`.

## Code Conventions

- **Indentation**: 4 spaces
- **Naming**: camelCase or snake_case for functions/variables; PascalCase for classes
- **Language**: British English (en-GB) — "colour", "behaviour", "modelling"
- When documenting a mizer S3 generic whose methods share a man page (combined with `@rdname`/`@name`), follow the steps in `.claude/skills/document-s3-generics.md`.
- When adding, moving or removing a species parameter default, follow `.claude/skills/species-param-defaults.md`. A default belongs to the rate setter that reads the parameter; only parameters that no single rate setter owns are defaulted centrally.
- When writing anything that integrates over the size grid — a summary or indicator function, a diagnostic derived from a rate, a rate setter with a size-dependent parameter — follow `.claude/skills/size-grid-integrals.md`. Each bin integral is performed in exactly one place, so a size-dependent factor is bin-averaged where its integral is performed and nowhere else; doing it twice is a silent uniform error.

## Testing

- Run `devtools::load_all()` before running tests
- Run only relevant tests with `devtools::test(filter = "pattern")`. Running all tests is too slow.
- **Never compare two rate arrays with `identical()`.** The freeze mechanism writes a `comment` attribute onto an array that has been set by hand, so two arrays holding the same numbers compare unequal whenever one of them has been through a setter. This reports a difference where there is none and has already produced a wrong conclusion about whether a setter argument works. Compare values: `all.equal()` after stripping attributes down to `dim`, or `all.equal(..., check.attributes = FALSE)`. The internal `different()` (`R/helpers.R:14`) is the in-package form of this.
- Use `expect_doppelganger()` (vdiffr) for plot tests, and snapshot tests for complex outputs
- When snapshot tests fail because values legitimately changed, run `testthat::snapshot_accept()` to promote the `.new.md` files into the canonical `.md` snapshots.
- Run `devtools::document()` after adding or changing exports
- Lint a file with `lintr::lint()`
- After editing C++ source: `devtools::clean_dll(); devtools::load_all()`
- After modifying the `MizerParams` or `MizerSim` class (new/removed slots, changes to `@rates_funcs`, etc.), follow the steps in `.claude/skills/upgrade-mizer-data.md`.
- Tests go in the file named after the R file that **defines** the function under test — not one named after the feature. Before adding a test, or any new test file, follow `.claude/skills/test-organisation.md`.
- Use the shared `_small` fixtures from `helper.R`; never a top-level `data()` call, as `testthat 3.x` shares `.GlobalEnv` across all tests. The fixtures give **3 species and 4 time steps**, so any hardcoded species or time index must be valid for that. For the full catalogue, for adding a fixture, and for the parallel and experimental-test settings, see `.claude/skills/test-fixtures.md`.

## Before Submitting

- After adding a new file under `R/`, add it to the `Collate:` field in `DESCRIPTION` (roxygen2 does not manage this automatically in this package).
- After adding or renaming exported functions, add them to the `appropriate` section in `pkgdown/_pkgdown.yml` so they appear on the reference page of the website.
- Update `NEWS.md` when adding features or fixing bugs.
- If the change can alter the behaviour of code users have already written — a renamed or removed argument, a changed default, a deprecation, a corrected bug that moves results — also record it in `inst/skills/upgrade-mizer-code/SKILL.md`, under the heading for the release being prepared, and add a row to its agent-only symptom index giving the symptom a user would actually see. Then regenerate (`build_cheatsheets()`), which updates `vignettes/upgrading.Rmd`. Do not edit that article directly.
- When creating a pull request, always include the summary of your changes in the PR body.
