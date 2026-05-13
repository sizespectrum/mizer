# get_initial_n ----
initial_n_params <- newMultispeciesParams(NS_species_params_gears, inter,
                                          info_level = 0)

test_that("get_initial_n returns steady-state initial abundances", {
    params <- initial_n_params
    n <- get_initial_n(params)
    expect_s3_class(n, "ArraySpeciesBySize")
    expect_equal(dim(n), dim(params@initial_n))
    expect_true(all(n >= 0))
    expect_equal(n, get_initial_n(params), ignore_attr = TRUE)
})

test_that("get_initial_n validates params and deprecates n0_mult", {
    expect_error(get_initial_n(1), "params argument must of type MizerParams")

    params <- initial_n_params
    lifecycle::expect_deprecated(n1 <- get_initial_n(params, n0_mult = 1),
                                 "deprecated")
    n2 <- get_initial_n(params)
    expect_equal(n1, n2, ignore_attr = TRUE)

    lifecycle::expect_deprecated(n3 <- get_initial_n(params, a = 1),
                                 "deprecated")
    expect_equal(n3, n2, ignore_attr = TRUE)
})
