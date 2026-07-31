This is a patch release (3.2.1) of the mizer package, which is already on
CRAN. It fixes how species and gear parameters are handled when they are
assigned: a size given as a weight can now actually be set on a model
specified by lengths, editing one of these parameter tables on its own no
longer validates it on every assignment, and `species_params<-()` no longer
errors on a list, matrix or object column. See NEWS.md for details.

## Test environments

### devtools::check_win_devel()

  TODO: run before submitting and record the result here.

### local Ubuntu, R 4.6.1

  0 errors ✔ | 0 warnings ✔ | 1 note ✖

  Compilation used the following non-portable flag(s):
     ‘-mno-omit-leaf-frame-pointer’

  This note is caused by a non-portable flag in the local R installation's
  compiler configuration, not by the package. It does not appear on the
  win-builder or mac-builder checks.

## Previous CRAN check results

The CRAN checks for the released version 3.2.0
(https://cran.rstudio.org/web/checks/check_results_mizer.html) are clean.

## Reverse dependencies

There are no reverse dependencies on CRAN.
