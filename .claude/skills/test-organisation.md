# Which file does a test belong in?

Use this skill when adding a test, moving one, or wondering whether a new test
file is warranted. The instinct to name a test file after the *feature* you are
working on is wrong here, and it is how the suite drifted out of shape before.

## Rule

A test goes in the test file named after the R file that **defines** the
function under test.

A test for `calc_selectivity()` — defined in `R/setFishing.R` — belongs in
`tests/testthat/test-setFishing.R`. Not `test-calc-selectivity.R`, not
`test-selectivity.R`. Find the home before writing the test:

```
grep -ln "^[ ]*<name> <- function" R/*.R
```

**Do not create a new test file just because a function or a feature is new.**
Create one only when a new file appears under `R/`. Conversely, every R file
containing logic should have a matching `test-<name>.R`.

The name must match exactly, including any suffix: `R/ArraySpeciesBySize-class.R`
is tested by `test-ArraySpeciesBySize-class.R`.

## The exception: genuinely cross-cutting tests

Some tests exercise the interplay of functions from several R files and so have
no single home. Those live in files named after the functionality they test.
This is a small, closed list — do not add to it without good reason:

- `test-analytic_results.R`, `test-analytic_transport.R` — rates and transport
  checked against closed-form solutions
- `test-single_species.R` — the whole pipeline on a one-species model
- `test-second_order_summary.R` — summary integrals under the `second_order_w` flag
- `test-backwards_compatibility.R` — legacy constructors and `R/deprecated.R`
- `test-reexports.R` — package-level re-exports

A test that mostly calls one function is not cross-cutting, even if it needs a
projection to set up. Ask which file you would edit if the test failed.

Some R files are deliberately without a test file of their own:
`deprecated.R` (covered by `test-backwards_compatibility.R`), `mizer-package.R`
(documentation only) and `RcppExports.R` (generated).

## Moving a test between files

Snapshots live in `_snaps/<test-file-name-without-the-test--prefix>.md`, keyed by
the `test_that()` description. Moving a test means moving its snapshot block to
the destination's `.md` file, or the snapshot is orphaned in one file and
missing in the other. vdiffr figures live in `_snaps/<name>/` directories the
same way.

If the destination file has top-level setup, check the moved test does not
depend on objects that only existed in its old file, and that its own
top-level objects do not collide with names already there.

## When an R file is split

Splitting an R file means splitting its test file the same way, so the
correspondence survives. Also add the new file to `Collate:` in `DESCRIPTION`
(roxygen2 does not manage this in mizer) and re-run `devtools::document()` — the
`% Please edit documentation in R/...` line of every moved function's `.Rd`
needs regenerating.
