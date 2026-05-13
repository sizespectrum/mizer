test_that("defaults_edition is a deprecated compatibility shim", {
    withr::local_options(mizer_defaults_edition = 1)

    lifecycle::expect_deprecated(result <- defaults_edition(), "deprecated")
    expect_identical(result, 2)
    expect_identical(getOption("mizer_defaults_edition"), 1)

    lifecycle::expect_deprecated(result <- defaults_edition(1), "deprecated")
    expect_identical(result, invisible(2))
})
