# get_initial_n ----
test_that("get_initial_n is working properly", {
    params <- newMultispeciesParams(NS_species_params_gears, inter, info_level = 0)
    n <- get_initial_n(params)
    no_sp <- nrow(params@species_params)
    for (i in 1:no_sp) {
        expect_true(all(n[i, params@w > params@species_params$w_max[i]] == 0))
        expect_true(all(n[i, params@w < params@species_params$w_min[i]] == 0))
    }
    # Check slope of all species is the same
    slopes <- rep(NA, no_sp)
    for (i in 1:no_sp) {
        n_idx <- which(n[i, ] != 0)
        slopes[i] <- (log(n[i, min(n_idx)]) - log(n[i, max(n_idx)])) / 
            (log(params@w[min(n_idx)]) - log(params@w[max(n_idx)]))
    }
    expect_equal(slopes, rep(slopes[1], no_sp))
    # Check that slopes = slope0
})
