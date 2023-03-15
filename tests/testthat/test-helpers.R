test_that("l2w works", {
    no_sp <- nrow(NS_species_params)
    # call with species_params
    expect_identical(l2w(2, NS_species_params), rep(0.08, no_sp))
    # call with params
    expect_identical(l2w(2, NS_params), rep(0.08, no_sp))
    # call with wrong 2nd argument
    expect_error(l2w(2, 4),
                 "The second argument must be either ")
    # call with wrong 1st argument
    expect_error(l2w("a", NS_species_params),
                 "l is not a numeric or integer vector")
    expect_error(l2w(1:2, NS_species_params),
                 "The length of 'l'")
})

test_that("w2l works", {
    no_sp <- nrow(NS_species_params)
    # call with species_params
    expect_identical(w2l(0.08, NS_species_params), rep(2, no_sp))
    # call with params
    expect_identical(w2l(0.08, NS_params), rep(2, no_sp))
    # call with wrong 1st argument
    expect_error(w2l("a", NS_species_params),
                 "w is not a numeric or integer vector")
    expect_error(w2l(1:2, NS_species_params),
                 "The length of 'w'")
    # w2l should do the inverse of l2w
    expect_equal(w2l(l2w(2, NS_species_params), NS_species_params),
                 rep(2, no_sp))
})
