# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(mizer)

# To keep CRAN check times down we run only a fast core subset of the test
# suite on CRAN, while the full suite still runs locally and on CI. The
# NOT_CRAN environment variable is set to "true" by devtools::test() and by
# R CMD check run locally, but is unset on CRAN's own machines (this is the
# same signal used by testthat::skip_on_cran()). The `filter` regex is matched
# against test file names with the "test-" prefix and ".R" suffix stripped.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    # Full suite locally and on CI
    test_check("mizer")
} else {
    # On CRAN: run only a fast core subset
    test_check("mizer",
               filter = "project$|MizerParams-class|newMultispeciesParams|rate_functions")
}
