test_that("l2w works", {
    no_sp <- nrow(NS_species_params)
    # call with species_params
    expect_identical(l2w(2, NS_species_params), rep(0.08, no_sp))
    # call with params
    expect_identical(l2w(2, NS_params), rep(0.08, no_sp))
    # call with wrong 2nd argument
    expect_error(l2w(2, 4),
                 "The second argument must be either ")
})

test_that("w2l works", {
    no_sp <- nrow(NS_species_params)
    # call with species_params
    expect_identical(w2l(0.08, NS_species_params), rep(2, no_sp))
    # call with params
    expect_identical(w2l(0.08, NS_params), rep(2, no_sp))
})
