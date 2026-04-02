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

test_that("constantEggRDI returns loss from the egg size bin", {
    expected <- {
        no_sp <- nrow(params@species_params)
        idx <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
        params@initial_n[idx] *
            (getEGrowth(params)[idx] + getMort(params)[idx] * params@dw[params@w_min_idx])
    }
    expect_equal(constantEggRDI(params, params@initial_n,
                                getEGrowth(params), getMort(params)),
                 expected)
})

test_that("mizerRDI integrates reproductive energy with erepro and egg size", {
    e_repro <- getERepro(params)
    expected <- 0.5 * drop((e_repro * params@initial_n) %*% params@dw) *
        params@species_params$erepro / params@w[params@w_min_idx]
    expect_equal(mizerRDI(params,
                          n = params@initial_n,
                          n_pp = params@initial_n_pp,
                          n_other = params@initial_n_other,
                          t = 0,
                          e_growth = getEGrowth(params),
                          mort = getMort(params),
                          e_repro = e_repro),
                 expected)
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
    sp$ricker_b <- seq_along(rdi) / 1000
    rdd <- RickerRDD(rdi, sp)
    expect_equal(rdd, rdi * exp(-sp$ricker_b * rdi))
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
