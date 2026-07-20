# Document an S3 generic with several methods

Use this skill when writing or editing roxygen2 documentation for a
mizer-defined **S3 generic whose methods share a single man page** (combined via
`@rdname`/`@name`). Following it avoids the most common failure here: the
"Documented arguments not in \usage" `R CMD check` error.

## Rule of thumb

Document **everything on the generic**. Leave each method with no roxygen beyond
`@rdname`, `@usage NULL` and `@export`.

## Where each argument is documented

- **Shared by all methods** → make it a formal of the generic and document it
  with an ordinary `@param`. Putting shared args in the generic's signature keeps
  `\usage` to a single short line (just the generic) and is what prevents the
  check error.
- **Used by only some methods** → do *not* give it a standalone `@param` (that
  re-introduces the check error). List it in a `\describe{}` block under
  `@param ...` on the generic.

Adding shared args to the generic's signature is safe: the generic body is just
`UseMethod()`, so its defaults are never evaluated, and `missing()` inside a
method still reflects the caller's actual call.

## Defaults

- **Class-dependent default** (e.g. `log_x` differs between size and time plots):
  give the arg *no* default in the generic signature, and explain the per-class
  default in its `@param` prose.
- **Object-dependent default** (e.g. `min_w = min(object@params@w)`): keep the
  arg under `@param ...` rather than in the generic signature.

## Avoid `@inheritParams`

Prefer inlining shared descriptions. On a multi-method page, `@inheritParams`
imports docs for the *union* of all method formals — including method-only ones —
as standalone `@param`, which re-triggers the check error.

## Exception: base-R generics

`plot`, `print`, `summary`, `as.data.frame` and other base-R generics cannot gain
new formals. Two strategies exist, depending on how much the methods' real
signatures vary. Prefer strategy 2 whenever the methods have more than a
couple of arguments each, or their signatures differ a lot (`plot` is the
example that forced this decision — see below); it's more verbose to set up
but has no CRAN-check gotchas and zero manual-sync burden afterwards.

### Strategy 1 (small/stable signatures: `print`, `summary`, `as.data.frame`, `str`)

Keep one shared page (`@name <generic>` on a `NULL` placeholder, or on one
method), give every method `@usage NULL`, and document arguments in a
`@param ...` `\describe{}` block. See the pitfalls below — this still requires
one manually-synced `\method{}` `@usage` line per method.

### Strategy 2 (large/varying signatures: `plot`)

Give every concrete method (`plot.ArraySpeciesBySize`, etc.) its own normal
roxygen doc block — real `@param` for exactly its own formals, `@keywords
internal` (hides it from the pkgdown reference index and softens
`checkDocFiles`/`checkRdContents` for that page), `@export`. Do **not** set
`@name`/`@rdname` or `@usage NULL` — let roxygen auto-derive `\usage` and
`\arguments` from the real function, the same as any other exported function
in the package. This is inherently correct and never goes stale, because nothing
is hand-transcribed.

Then add a *separate*, purely conceptual page (`@name plot` on a `NULL`
placeholder) with **no `@param` at all** — describe arguments in prose inside
`@details` (e.g. a `\describe{}` list), not via `@param`. A page with neither
`\usage` nor `\arguments` is never inspected by `codoc()`/`checkDocFiles()` —
confirmed by several pre-existing mizer pages with this exact shape (e.g.
`plotting_functions.Rd`, `mizer-package.Rd`) that have never been flagged.
Link to the individual method pages with `[plot.ArraySpeciesBySize()]` etc. so
users can find exact per-class signatures.

This is what backs the current `plot.Rd` (conceptual) plus
`plot.ArraySpeciesBySize.Rd`, `plot.ArrayTimeBySpecies.Rd`,
`plot.ArrayTimeBySpeciesBySize.Rd`, `plot.ArrayResourceBySize.Rd`,
`plot.ArrayTimeByResourceBySize.Rd` (individual, internal) — use that as the
template for adding another varying-signature method later.

### Strategy 1 pitfalls (still relevant for `print`/`summary`/`str`/`as.data.frame`)

All of their method arguments go in the `@param ...` `\describe{}`
block on the page, and the methods use `@usage NULL`.

Because every method on the page uses `@usage NULL` (and, for `print`/`summary`/
`as.data.frame`/`str`, the shared doc block is attached to a `NULL` placeholder
rather than a real function), roxygen2 never auto-generates a `\usage{}` section
for these pages at all. An Rd file with `\arguments` but no `\usage` triggers the
`R CMD check` NOTE "Rd files without \usage".

**Do not "fix" this by writing a `@usage` line for the bare generic** (e.g.
`@usage print(x, ...)`). Mizer doesn't define a function literally named
`print`/`plot`/`summary`/`as.data.frame`/`str` — only methods like
`print.ArraySpeciesBySize`. `tools::codoc()` (the "code/documentation
mismatches" check `R CMD check` runs) looks up an object of that exact name and,
finding none, raises a WARNING — worse than the NOTE it replaced — "Functions or
methods with usage in Rd file '<name>.Rd' but not in code: '<name>'".

The correct fix is one `\method{generic}{class}(...)` line per aliased method,
reproducing that method's *real, full* formal list (names, order and defaults)
exactly as defined in `R/`, e.g.:

```r
#' @usage
#' \method{print}{ArraySpeciesBySize}(x, ...)
#' \method{print}{ArrayTimeBySpecies}(x, ...)
#' \method{str}{MizerSim}(object, max.level = NA, ...)
```

`tools::codoc()` resolves `\method{generic}{class}` via `getS3method(generic,
class)` and diffs the formal list against it — so every alias on the page needs
its own line, kept in sync whenever that method's signature changes. Verify
with `tools::checkRd("man/<name>.Rd")` (catches malformed Rd) *and* a real
`tools::codoc()` run against an **installed** copy of the package (catches
signature mismatches; `checkRd` alone won't):

```r
devtools::document()
install.packages(".", repos = NULL, type = "source", lib = "/tmp/rlib")
tools::codoc("mizer", lib.loc = "/tmp/rlib")
```

A clean result is `length(res) == 0` with `attr(res,
"functions_in_usages_not_in_code")` empty.

**This full-signature `\usage` breaks the "nest under `...`" rule above for any
argument that is a real formal of at least one method.** The "Where each
argument is documented" rule (put method-only args in a `\describe{}` under
`@param ...`) only works when `\usage` hides individual method signatures
(`@usage NULL` on every method, one short generic-only line at the top). Once
`\usage` spells out each method's real arguments (as it must here, for
`codoc()`), `tools::checkDocFiles()` (run by `devtools::check_man()`, and by
`R CMD check` as "Undocumented arguments in documentation object") requires a
matching **top-level** `@param <name>` for every argument name that appears in
*any* `\usage` line — a name buried in a nested `\describe{}` doesn't count.
(This is exactly why `plot` moved to strategy 2: with 5 methods and up to 16
varying arguments, flattening every one to a top-level `@param` on a single
shared page became unreadable — see strategy 2 above.) For pages small enough
to stay on strategy 1: give every real argument (not just ones shared by every
method) its own flat `@param`, and say which methods it applies to in the
prose (e.g. "For `ArrayTimeBySpecies` methods, ..."), rather than nesting it.
Keep `@param ...` itself only for genuine passthrough (not a named formal of
any method). Verify with `devtools::check_man()` — it should report nothing.

## The `param_object_dots` template

Most rate-function generics use `@template param_object_dots`. That template
already supplies `@param object` and the full `@param ...` `\describe{}` block
covering `n`, `n_pp`, `n_other`, `t` (MizerParams) and `time_range`, `drop`
(MizerSim). You only need to add `@param` entries for *extra* shared args (like
`power`) that are not in the template.

## After editing

Run `devtools::document()` to regenerate `man/` from the roxygen2 source, then
check the updated `.Rd` to confirm:
- `\usage{}` contains every `@param`-documented argument.
- Both methods appear as `\alias{}` entries on the shared page.
