
test_that("getRequiredRDD works for single species model", {
    params <- single_sp_params
    # In a steady state model, getRequiredRDD should return the same as getRDD
    # because the reproduction matches the required amount to maintain the steady state.
    rdd <- getRequiredRDD(params)
    rdd_actual <- getRDD(params)

    expect_equal(rdd, rdd_actual)
})

test_that("getRequiredRDD works for community model", {
    params <- community_params_small
    # In community model, reproduction is constant
    rdd <- getRequiredRDD(params)
    rdd_actual <- getRDD(params)

    # Check if they are close enough
    expect_equal(rdd, rdd_actual)
})

test_that("getRequiredRDD handles diffusion", {
    # Create a model with diffusion
    params <- single_sp_params

    # Add diffusion
    ext_diff <- params@ext_diffusion
    ext_diff[] <- 1e9 * params@w
    params <- setExtDiffusion(params, ext_diffusion = ext_diff)

    # Update initial_n to be the steady state solution with this diffusion
    # We need to recalculate it using get_steady_state_n

    # Instead of fully replicating newSingleSpeciesParams logic, let's use the fact that
    # getRequiredRDD should make the current state a steady state *at the boundary*.

    # If we simply update reproduction to match getRequiredRDD, then the boundary flux matches reproduction.

    # So let's update erepro
    rdd_req <- getRequiredRDD(params)
    params@species_params$erepro <- params@species_params$erepro * rdd_req / getRDI(params)

    # Now getRDD(params) should match rdd_req
    expect_equal(getRDD(params), rdd_req)

    # Let's stick to the test that it does not error and returns a numeric vector of correct length
    expect_length(rdd_req, 1)
    expect_type(rdd_req, "double")
    expect_true(rdd_req > 0)

    # verifying that it is different from the no-diffusion case
    params_no_diff <- single_sp_params
    expect_true(rdd_req != getRequiredRDD(params_no_diff))
})

test_that("getRequiredRDD uses total diffusion", {
    params <- single_sp_params
    species <- params@species_params$species[1]
    params@use_predation_diffusion <- TRUE

    idx <- params@w_min_idx[species]
    next_idx <- idx + 1
    dw <- params@dw[idx]
    diffusion <- getDiffusion(params)

    expected <- initialN(params)[species, idx] *
        (getEGrowth(params)[species, idx] +
             getMort(params)[species, idx] * dw) +
        0.5 * (diffusion[species, idx] * initialN(params)[species, idx] -
                   diffusion[species, next_idx] *
                   initialN(params)[species, next_idx]) / dw

    expect_equal(getRequiredRDD(params)[species], expected,
                 ignore_attr = TRUE)
})

test_that("getRequiredRDD keeps egg density constant with diffusion", {
    params <- single_sp_params
    species <- params@species_params$species[1]
    n <- params@species_params[species, "n"]
    ext_diffusion(params)[species, ] <- 0.1 * params@w^(n + 1)
    params <- steadySingleSpecies(params, species = species)
    params@species_params$constant_reproduction <- getRequiredRDD(params)
    params <- setRateFunction(params, "RDD", "constantRDD")

    sim <- project(params, t_max = 1, t_save = 1, progress_bar = FALSE)
    idx <- params@w_min_idx[species]
    expect_equal(finalN(sim)[species, idx], initialN(params)[species, idx],
                 tolerance = 1e-12)
})

test_that("getRequiredRDD works with predictor-corrector project method", {
    params <- single_sp_params
    species <- params@species_params$species[1]
    n <- params@species_params[species, "n"]
    ext_diffusion(params)[species, ] <- 0.1 * params@w^(n + 1)
    params <- steadySingleSpecies(params, species = species)
    params@species_params$constant_reproduction <- getRequiredRDD(params)
    params <- setRateFunction(params, "RDD", "constantRDD")

    sim <- project(params, t_max = 1, t_save = 1,
                   method = "predictor-corrector", progress_bar = FALSE)
    idx <- params@w_min_idx[species]
    expect_equal(finalN(sim)[species, idx], initialN(params)[species, idx],
                 tolerance = 1e-12)
})

test_that("required recruitment follows a complete rate replacement", {
    rlang::local_bindings(
        required_full_rates = function(params, ...) {
            r <- mizerRates(params, ...)
            r$mort <- 2 * r$mort
            r
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    baseline <- getRequiredRDD(p)
    egg <- cbind(seq_len(nrow(p@initial_n)), p@w_min_idx)
    extra <- getMort(p)[egg] * p@initial_n[egg] * p@dw[p@w_min_idx]
    p <- setRateFunction(p, "Rates", "required_full_rates")
    expect_equal(unname(getRequiredRDD(p)), unname(baseline + extra))

    fns <- projectRateFunctions(p)
    r <- fns$Rates(p, n = p@initial_n, n_pp = p@initial_n_pp,
                   n_other = p@initial_n_other, effort = p@initial_effort,
                   rates_fns = fns)
    expect_equal(getRequiredRDD(p, rates = r), getRequiredRDD(p))
})

test_that("required recruitment follows full-rate S3 dispatch", {
    method <- function(params, ...) {
        r <- NextMethod()
        r$mort <- 2 * r$mort
        r
    }
    method_name <- paste0("projectRates.", "required_test")
    assign(method_name, method, envir = .GlobalEnv)
    withr::defer(rm(list = method_name, envir = .GlobalEnv))
    p <- NS_params_small
    baseline <- getRequiredRDD(p)
    egg <- cbind(seq_len(nrow(p@initial_n)), p@w_min_idx)
    extra <- getMort(p)[egg] * p@initial_n[egg] * p@dw[p@w_min_idx]
    class(p) <- c("required_test", class(p))
    expect_equal(unname(getRequiredRDD(p)), unname(baseline + extra))
})

test_that("egg balance does not recursively request ordinary RDI", {
    rlang::local_bindings(
        required_rdi = function(params, ...) getRequiredRDD(params),
        .env = .GlobalEnv
    )
    p <- NS_params_small
    baseline <- getRequiredRDD(p)
    p <- setRateFunction(p, "RDI", "required_rdi")
    expect_equal(getRequiredRDD(p), baseline)
    expect_equal(as.numeric(getRDI(p)), as.numeric(baseline))
})

test_that("circular complete-rate dependencies give an actionable error", {
    rlang::local_bindings(
        recursive_required_rates = function(params, ...) {
            mizerRates(params, ...)
        },
        recursive_required_rdi = function(params, ...) getRequiredRDD(params),
        .env = .GlobalEnv
    )
    p <- setRateFunction(NS_params_small, "Rates", "recursive_required_rates")
    p <- setRateFunction(p, "RDI", "recursive_required_rdi")
    expect_error(getRequiredRDD(p), "Recursive full-rate evaluation")
})
