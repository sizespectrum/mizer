test_that("defaults_edition gets and sets the active edition", {
    old <- getOption("mizer_defaults_edition")
    on.exit(options(mizer_defaults_edition = old), add = TRUE)
    options(mizer_defaults_edition = NULL)

    expect_message(result <- defaults_edition(200),
                   "Mizer parameter defaults are now at edition 200")
    expect_identical(result, invisible(200))
    expect_identical(defaults_edition(), 200)
})

test_that("defaults_edition validates the supplied edition", {
    expect_error(defaults_edition("2"), "not a numeric or integer vector")
    expect_error(defaults_edition(0), "not greater than or equal to 1")
})
