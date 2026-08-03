# Test fixtures and keeping the suite fast

Use this skill when writing a test that needs a model to work with, when adding
a shared fixture to `helper.R`, or when a test turns out to be slow.

## Use the shared `_small` fixtures

`helper.R` defines small objects used across all test files. They are named with
the `_small` suffix so they cannot shadow the built-in package datasets.

- **`NS_params_small`** — 3 species (Sprat=1, Herring=2, Cod=3), 20 size bins,
  3 gears (Industrial effort=0, Pelagic effort=1, Otter effort=0.5). Sprat is
  always unfished (Industrial gear, zero effort). Used by ~60 test files.
- **`NS_sim_small`** — `project(NS_params_small, t_max = 3, t_save = 1)`,
  giving 4 time steps at t=0,1,2,3.
- **`NS_species_params_small`**, **`NS_species_params_gears_small`**,
  **`inter_small`** — 3-species subsets of the package datasets
  (rows 1, 4, 11 = Sprat, Herring, Cod).
- **`NS_params_default_small`** — the same 3 species but with the default
  `no_w`, for tests that must not depend on the 20-bin grid.
- **`NS_params_cod_small`** — single-species (Cod) model.
- **`single_sp_params`**, **`trait_params_small`**, **`trait_params_2sp`**,
  **`community_params_small`** — the wrapper-constructor models.
- **`example_params()`** — a function building a 3-species model with differing
  egg sizes, length-weight parameters, a multi-gear `gear_params` and diffusion
  on one species.

**Fixtures give 3 species and 4 time steps.** Any hardcoded species index
(`[12, ]`) or time index (`times[10]`) must be valid for that.

Never use a top-level `data()` call in a test file: `testthat 3.x` shares
`.GlobalEnv` across all tests, so it leaks into every other file.

## Adding a fixture

Build it from scratch. Do **not** load a fixture from a saved `.rds`: a
serialized `MizerParams` is a snapshot of the class as it was when saved, so it
would silently stop reflecting changes to a slot, to a rate setter's default or
to a constructor — and the failure mode is tests that keep *passing* against an
obsolete object. That is the problem `upgradeParams()` exists to solve for saved
objects, and it is not one worth importing into the fixtures. Building from
scratch also means every run exercises the constructors for free.

Wrap it in `delayedAssign()` unless nearly every file needs it:

```r
delayedAssign("trait_params_small", suppressMessages(newTraitParams()))
```

Every parallel worker is a separate process and sources `helper.R` in full, so
an eager fixture is built once per worker whether or not that worker's files use
it. `delayedAssign()` is transparent — test files still refer to the bare name —
and a worker then only builds what it actually touches. Keep the
`suppressMessages()` inside the promise, so constructor output cannot leak into
whichever test happens to force it first.

## Running the suite

`devtools::load_all()` first, then `devtools::test(filter = "pattern")`. Running
everything is rarely what you want while iterating.

The suite runs in parallel (`Config/testthat/parallel: true`). testthat uses
**2 workers by default whatever the machine has** — that default is deliberate,
since CRAN caps checks at 2 cores. Raise it locally with `TESTTHAT_CPUS=8` or
`options(Ncpus = 8)`, though the gain tails off quickly once the longest single
file sets the floor.

`Config/testthat/start-first` in `DESCRIPTION` lists the longest files so they
start immediately instead of becoming a lone tail at the end of the run. If a
file becomes slow, add it there.

## Slow tests of experimental features

Tests of functions marked `lifecycle::badge("experimental")` may be put behind an
opt-in flag when they are slow. This currently covers `test-steadyNewton.R` and
`test-getLimitCycleSim.R`, which together were a third of the suite's runtime.

Gate each block with the helper from `helper.R`:

```r
test_that("...", {
    skip_unless_experimental()
    ...
})
```

and keep expensive top-level setup in a `delayedAssign()`, or it still runs when
every test in the file skips.

Run them explicitly with:

```
MIZER_TEST_EXPERIMENTAL=true Rscript -e 'devtools::test(filter = "steadyNewton")'
```

The `check-standard` and `basic_check` workflows set the variable, so a full
R CMD check still covers this code. **Remove the gate when the function leaves
experimental status** — otherwise the tests quietly stop running for good.
