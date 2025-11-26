# choose an example params object
params <- NS_params
sp <- params@species_params
rdi <- getRDI(params)

test_that("constantEggRDI() keeps egg density constant", {
    # We set the reproduction rate functions
    params <- setRateFunction(params, "RDI", "constantEggRDI")
    params <- setRateFunction(params, "RDD", "noRDD")
    # Now the egg density stays fixed no matter how we fish
    sim <- project(params, t_max = 1, effort = 1)
    # Check that indeed the egg densities have not changed
    no_sp <- nrow(params@species_params) # number of species
    # Hacky shortcut to access the correct element of a 2D array using 1D notation
    idx <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
    expect_equal(finalN(sim)[idx], initialN(params)[idx])
})

test_that("BevertonHoltRDD works", {
    rdd <- BevertonHoltRDD(rdi, sp)
    expect_identical(rdd, rdi / (1 + rdi/sp$R_max))
})

test_that("BevertonHoltRDD checks for R_max column", {
    sp_no_rmax <- sp
    sp_no_rmax$R_max <- NULL
    expect_error(BevertonHoltRDD(rdi, sp_no_rmax),
                 "The R_max column is missing in species_params.")
})

test_that("RickerRDD works", {
    expect_error(RickerRDD(rdi, sp),
                 "The ricker_b column is missing in species_params")
    sp$ricker_b <- 0
    rdd <- RickerRDD(rdi, sp)
    expect_identical(rdd, rdi)
})

test_that("SheperdRDD works", {
    expect_error(SheperdRDD(rdi, sp),
                 "The species_params dataframe must contain columns sheperd_b and sheperd_c.")
    sp$sheperd_b <- 1/sp$R_max
    sp$sheperd_c <- 1
    rdd <- SheperdRDD(rdi, sp)
    expect_equal(rdd, BevertonHoltRDD(rdi, sp))
})

test_that("constantRDD returns constant_reproduction values", {
    sp$constant_reproduction <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200)
    rdd <- constantRDD(rdi, sp)
    expect_identical(rdd, sp$constant_reproduction)
})

test_that("noRDD returns rdi unchanged", {
    rdd <- noRDD(rdi, sp)
    expect_identical(rdd, rdi)
})

test_that("w_repro_max is rounded down to w grid for new versions", {
    # Create a new params object with current version
    params_new <- NS_params
    # Set a w_repro_max that doesn't fall exactly on a grid point
    # Use a value between two grid points
    w_grid <- params_new@w
    mid_idx <- length(w_grid) %/% 2
    original_w_repro_max <- (w_grid[mid_idx] + w_grid[mid_idx + 1]) / 2
    params_new@species_params$w_repro_max <- rep(original_w_repro_max, 
                                                  nrow(params_new@species_params))
    
    # Call setReproduction which should round down w_repro_max
    params_new <- setReproduction(params_new)
    
    # Check that w_repro_max was rounded down to a grid point
    for (i in seq_len(nrow(params_new@species_params))) {
        expect_true(params_new@species_params$w_repro_max[i] %in% params_new@w)
        expect_true(params_new@species_params$w_repro_max[i] <= original_w_repro_max)
    }
    
    # Verify it's the largest grid point <= original value
    expected_w_repro_max <- w_grid[mid_idx]
    expect_equal(params_new@species_params$w_repro_max[1], expected_w_repro_max)
})

test_that("w_repro_max rounding only applies to new versions", {
    # Create an old params object by setting mizer_version to an old value
    params_old <- NS_params
    params_old@mizer_version <- as.package_version("2.5.3.9000")
    
    # Set a w_repro_max that doesn't fall exactly on a grid point
    w_grid <- params_old@w
    mid_idx <- length(w_grid) %/% 2
    original_w_repro_max <- (w_grid[mid_idx] + w_grid[mid_idx + 1]) / 2
    params_old@species_params$w_repro_max <- rep(original_w_repro_max, 
                                                  nrow(params_old@species_params))
    
    # Call setReproduction - should NOT round for old versions
    params_old <- setReproduction(params_old)
    
    # Check that w_repro_max was NOT rounded and remains as set
    # For old versions, it should keep the off-grid value
    expect_equal(params_old@species_params$w_repro_max[1], original_w_repro_max)
})
