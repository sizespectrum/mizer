This is a minor release (3.2.0) of the mizer package, which is already on
CRAN. It overhauls how species and resource parameters are set, makes the
extension framework composable regardless of load order, adds a new
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

The current CRAN checks for the released version 3.1.0
(https://cran.rstudio.org/web/checks/check_results_mizer.html) report a NOTE on
the r-devel Debian flavors about several .Rd files (`as.data.frame.Rd`,
`plot.Rd`, `print.Rd`, `str.Rd`, `summary.Rd`) with documented arguments but no
`\usage` section. Re-checking the current version on a Debian r-devel platform
with rhub (the `gcc14` container) gives `Status: OK` with
`checking Rd \usage sections ... OK`, so this NOTE no longer occurs.

## Reverse dependencies

There are no reverse dependencies on CRAN.
