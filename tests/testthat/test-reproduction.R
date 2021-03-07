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