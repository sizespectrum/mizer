This is a resubmission. The win-builder incoming-pretest log for the previous
submission attempt reported a NOTE (`Rd files without \usage:
'as.data.frame.Rd' 'plot.Rd' 'print.Rd' 'str.Rd' 'summary.Rd'`) on the Debian
r-devel/gcc-16 platform. This has now been fixed directly by adding explicit
`@usage` tags to these five shared man pages; see "Previous CRAN check
results" below for details.

This is otherwise a minor release (3.2.0) of the mizer package, which is
already on CRAN. It overhauls how species and resource parameters are set,
makes the extension framework composable regardless of load order, adds a new
`adjustSizeGrid()` function, and includes a range of smaller improvements and
bug fixes. See NEWS.md for details.

## Test environments

### devtools::check_win_devel()

  Status: OK

### local Ubuntu, R 4.6.1

  0 errors ✔ | 0 warnings ✔ | 1 note ✖

  Compilation used the following non-portable flag(s):
     ‘-mno-omit-leaf-frame-pointer’

  This note is caused by a non-portable flag in the local R installation's
  compiler configuration, not by the package. It does not appear on the
  win-builder or mac-builder checks above.

## Previous CRAN check results

This same NOTE was already visible in the CRAN checks for the released version
3.1.0 (https://cran.rstudio.org/web/checks/check_results_mizer.html), but only
on the r-devel Debian flavors there too — it does not show up in
`devtools::check_win_devel()` (Windows only) or on rhub's Debian `gcc14`
container.

The root cause has now been fixed directly: these five pages document methods
of base-R generics (`plot`, `print`, `summary`, `as.data.frame`, `str`), and
every method on each page uses `@usage NULL` (required since base generics
can't gain new formals), which left roxygen2 generating no `\usage{}` section
at all. Each page now carries an explicit `@usage` tag with the generic's own
signature, e.g. `print(x, ...)`. Verified locally with `tools::checkRd()` and
`devtools::check_man()`; since neither `check_win_devel()` nor a local Ubuntu
check can reproduce the original platform/compiler combination, this fix could
not be re-confirmed on that exact combination before submission.

## Reverse dependencies

There are no reverse dependencies on CRAN.
