This is a resubmission. It fixes the NOTE and the subsequent WARNING raised
against the previous submission attempts, both about `\usage`/`\arguments` in
`as.data.frame.Rd`, `plot.Rd`, `print.Rd`, `str.Rd` and `summary.Rd`.

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

The original NOTE was already visible in the CRAN checks for the
released version 3.1.0
(https://cran.rstudio.org/web/checks/check_results_mizer.html), on the
r-devel Debian flavors and this has now been fixed.

## Reverse dependencies

There are no reverse dependencies on CRAN.
