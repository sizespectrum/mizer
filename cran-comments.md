## Test environments

### devtools::check_win_devel()

  Status: 2 NOTEs
  
devtools::check_win_devel()

### local Ubuntu, R 4.6.1

  0 errors ✔ | 0 warnings ✔ | 3 notes ✖
 
As above with one extra note:

Compilation used the following non-portable flag(s):
     ‘-mno-omit-leaf-frame-pointer’

This note is caused by a non-portable flag in the local R installation's
compiler configuration, not by the package. It does not appear on the
win-builder or mac-builder checks.

## Previous CRAN check results

The CRAN checks for the released version 3.3.0
(https://cran.rstudio.org/web/checks/check_results_mizer.html) are clean.

## Reverse dependencies

There are no reverse dependencies on CRAN.
