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

test_that("a model records the edition active when it was created", {
    params <- withr::with_options(
        list(mizer_defaults_edition = 2),
        newMultispeciesParams(NS_species_params_small, no_w = 20,
                              info_level = 0))
    expect_identical(defaults_edition(params), 2)
    # It keeps it when another edition becomes active
    withr::local_options(mizer_defaults_edition = 1)
    expect_identical(defaults_edition(params), 2)
    expect_identical(defaults_edition(NS_params_small), 1)
    # A simulation has the edition of its model
    sim <- project(params, t_max = 0.1, t_save = 0.1)
    expect_identical(defaults_edition(sim), 2)
})

test_that("a model that has not recorded an edition is of edition 1", {
    params <- NS_params_small
    params$defaults_edition <- NULL
    expect_identical(defaults_edition(params), 1)
})

test_that("saved S4 models without an edition are of edition 1", {
    sim <- readRDS("assets/simc.1.0.rds")
    params <- methods::slot(sim, "params")
    expect_true(isS4(sim))
    expect_true(isS4(params))
    expect_false(methods::.hasSlot(params, "defaults_edition"))

    # Read the legacy objects directly, before any validation or upgrade.
    withr::local_options(mizer_defaults_edition = 2)
    expect_identical(defaults_edition(params), 1)
    expect_identical(defaults_edition(sim), 1)
    expect_identical(defaults_edition(), 2)
})

test_that("defaults_edition<- changes the recorded edition and nothing else", {
    params <- NS_params_small
    defaults_edition(params) <- 2
    expect_identical(defaults_edition(params), 2)
    changed <- names(params)[!mapply(identical, params, NS_params_small)]
    expect_setequal(changed, c("defaults_edition", "time_modified"))
    expect_error(defaults_edition(params) <- 0,
                 "not greater than or equal to 1")
    expect_error(defaults_edition(params) <- c(1, 2), "not a number")
})

test_that("with_defaults_edition restores the active edition", {
    withr::local_options(mizer_defaults_edition = 1)
    expect_identical(with_defaults_edition(2, defaults_edition()), 2)
    expect_identical(defaults_edition(), 1)
    expect_error(with_defaults_edition(2, stop("oops")), "oops")
    expect_identical(defaults_edition(), 1)
})
