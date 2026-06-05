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
new formals. So all of their method arguments go in the `@param ...` `\describe{}`
block on the page, and the methods use `@usage NULL`.
