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

### Documenting S3 generics with several methods

When several methods of a mizer-defined S3 generic share a man page (combined
with `@rdname`/`@name`), document everything on the **generic** and leave the
methods with no roxygen beyond `@rdname`, `@usage NULL` and `@export`:

- The generic's **function signature owns all arguments shared by all methods**.
  Document those with ordinary `@param`. This keeps `\usage` short (one generic
  signature) and avoids the "Documented arguments not in \usage" `R CMD check`
  error.
- Arguments used by **only some** methods are listed in a `\describe{}` block
  under `@param ...` on the generic (not as standalone `@param`, which would
  re-introduce the check error).
- Adding shared args to the generic is safe: the generic body is just
  `UseMethod()`, so its defaults are never evaluated and `missing()` in a method
  still reflects the caller's call.
- Where a shared arg's default is **class-dependent** (e.g. `log_x` differs
  between size and time plots), give it no default in the generic signature and
  explain the per-class default in its `@param` prose. Where a default is
  **object-dependent** (e.g. `min_w = min(object@params@w)`), keep that arg
  under `@param ...` instead.
- Prefer inlining shared descriptions over `@inheritParams`: on a multi-method
  page `@inheritParams` imports docs for the *union* of all method formals
  (including method-only ones) as standalone `@param`, which re-triggers the
  check error.
- **Exception:** base-R generics (`plot`, `print`, `summary`, `as.data.frame`)
  cannot gain formals, so all their method arguments go in the `@param ...`
  `\describe{}` block. For the `@usage`, reference the **generic**
  (`@usage plot(x, ...)`), not a `\method{generic}{class}(...)` line, and keep
  every method at `@usage NULL`. A fabricated `\method{}` usage that omits the
  method's real formals triggers an `R CMD check` *codoc* warning ("Codoc
  mismatches"), because codoc compares each `\method{}` usage against the actual
  method signature. Referencing the base generic sidesteps this: codoc skips it
  because the generic is not a package object. (A `\method{}` usage is only safe
  when it lists the method's real formals, e.g. `print` methods that genuinely
  take just `(x, ...)`.)

## Testing

- Use `expect_doppelganger()` (vdiffr) for plot tests
- Use snapshot tests for complex outputs
- Run `devtools::document()` after adding or changing exports
- Run `devtools::load_all()` before running tests
- After modifying the `MizerParams` or `MizerSim` class (new/removed slots, changes to `@rates_funcs`, etc.), follow the steps in `.claude/skills/upgrade-mizer-data.md`. Claude Code users can invoke this as `/upgrade-mizer-data`.

## Before Submitting

- After adding a new file under `R/`, add it to the `Collate:` field in `DESCRIPTION` (roxygen2 does not manage this automatically in this package).
- After adding or renaming exported functions, add them to the `appropriate` section in `pkgdown/_pkgdown.yml` so they appear on the reference page of the website.
- Update `NEWS.md` when adding features or fixing bugs.
