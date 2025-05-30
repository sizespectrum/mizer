test_that("age_mat_vB works with MizerParams object", {
    params <- NS_params
    ages <- age_mat_vB(params)

    expect_named(ages, params@species_params$species)
    expect_type(ages, "double")
    expect_true(all(ages > 0, na.rm = TRUE))
})

test_that("age_mat_vB works with species_params dataframe", {
    sp <- NS_params@species_params
    ages <- age_mat_vB(sp)

    expect_named(ages, sp$species)
    expect_type(ages, "double")
    expect_true(all(ages > 0, na.rm = TRUE))
})

test_that("age_mat_vB handles missing parameters correctly", {
    sp <- NS_params@species_params
    # Test missing k_vb
    sp_no_k <- sp
    sp_no_k$k_vb <- NULL
    ages <- age_mat_vB(sp_no_k)
    expect_true(all(is.na(ages)))

    # Test missing w_inf but with w_max
    sp_no_winf <- sp
    sp_no_winf$w_inf <- NULL
    ages <- age_mat_vB(sp_no_winf)
    expect_false(any(is.na(ages)))

    # Test missing b parameter
    sp_no_b <- sp
    sp_no_b$b <- NULL
    ages <- age_mat_vB(sp_no_b)
    expect_false(any(is.na(ages)))
    expect_true(all(ages > 0))
})

test_that("age_mat_vB throws error for invalid input", {
    expect_error(
        age_mat_vB("not a valid input"),
        "The first argument must be either a MizerParams object or a species_params data frame."
    )
})

test_that("age_mat works with MizerParams object", {
    params <- NS_params
    ages <- age_mat(params)

    expect_named(ages, params@species_params$species)
    expect_type(ages, "double")
    expect_true(all(ages > 0))
    expect_true(all(is.finite(ages)))
})

test_that("age_mat throws error for invalid input", {
    expect_error(age_mat("not a valid input"))
})
