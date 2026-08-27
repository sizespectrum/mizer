## Test environments

### devtools::check_win_release()

  Status: OK

### local Ubuntu, R 4.6.1

  0 errors ✔ | 0 warnings ✔ | 1 note ✖

  Compilation used the following non-portable flag(s):
     ‘-mno-omit-leaf-frame-pointer’

  This note is caused by a non-portable flag in the local R installation's
  compiler configuration, not by the package. It does not appear on the
  win-builder or mac-builder checks.

## The call to `attach()`

`R CMD check --as-cran` reports

```
Found the following calls to attach():
File 'mizer/R/registerExtensions.R':
  attach(NULL, name = .extensionClassEnvName, warn.conflicts = FALSE)
```

This one is deliberate. Mizer builds a chain of S4 marker classes at run time,
one per loaded extension package, so that extensions compose in load order.
That metadata cannot live in a package namespace, which is locked and so cannot
support the removal and re-parenting that rebuilding the chain needs. Nor can
it live in `.GlobalEnv`, which `cleanEx()` empties before every example: that
destroyed the classes and failed the examples of every extension package. An
environment that mizer merely holds privately does not serve either, because
`methods` resolves a class whose package slot names an extension package by an
inheriting lookup that runs through `.GlobalEnv` and the search path but never
through mizer's own namespace.

The attached environment is named `mizer:extension-classes`, holds nothing but
S4 class metadata, and masks nothing. It is attached only when a dispatching
extension package is loaded, so a session using base mizer never acquires it.

## Previous CRAN check results

The CRAN checks for the released version 3.2.1
(https://cran.rstudio.org/web/checks/check_results_mizer.html) are clean.

## Reverse dependencies

There are no reverse dependencies on CRAN.
