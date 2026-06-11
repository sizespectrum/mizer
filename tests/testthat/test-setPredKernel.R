params <- NS_params_small
no_sp <- nrow(params@species_params)

## setPredKernel ----
test_that("setPredKernel works", {
    expect_unchanged(setPredKernel(params), params)
    expect_unchanged(setPredKernel(params, pred_kernel = NULL), 
                     params)
    params@species_params$pred_kernel_type <- "box"
    params@species_params$ppmr_min <- 2
    expect_error(setPredKernel(params), 
                 "missing from the parameter dataframe: ppmr_max")
    params@species_params$ppmr_max <- 4
    p2 <- setPredKernel(params)
    pred_kernel <- 2 * getPredKernel(params)
    expect_error(setPredKernel(params, pred_kernel[1:2, ]),
                 "incorrect number of dimensions")
    expect_error(setPredKernel(params, pred_kernel - 1),
                 "pred_kernel >= 0 are not true")
    p2 <- setPredKernel(params, pred_kernel)
    expect_equal(p2@pred_kernel, pred_kernel, ignore_attr = TRUE)
    expect_identical(p2@pred_kernel, getPredKernel(p2))
})

test_that("Comment works on pred_kernel", {
    params <- NS_params_small
    # if no comment, it is set automatically
    pred_kernel <- getPredKernel(params)
    params <- setPredKernel(params, pred_kernel = pred_kernel)
    expect_identical(comment(params@pred_kernel), "set manually")
    
    # comment is stored
    comment(pred_kernel) <- "test"
    params <- setPredKernel(params, pred_kernel = pred_kernel)
    expect_identical(comment(params@pred_kernel), "test")
    
    # if no comment, previous comment is kept
    comment(pred_kernel) <- NULL
    params <- setPredKernel(params, pred_kernel = pred_kernel)
    expect_identical(comment(params@pred_kernel), "test")
    
    # no message when nothing changes
    expect_message(setPredKernel(params), NA)
    # but message when a change is not stored due to comment
    beta <- params@species_params$beta
    params@species_params$beta <- 1
    expect_message(setPredKernel(params),
                   "You have set a custom predation kernel")
    # Can reset
    params@species_params$beta <- beta
    p <- setPredKernel(params, reset = TRUE)
    expect_equal(p@pred_kernel, pred_kernel)
    expect_warning(setPredKernel(params, pred_kernel = pred_kernel,
                                    reset = TRUE),
                   "Because you set `reset = TRUE`, the")
})

# getPredKernel ----
test_that("getPredKernel has correct dimnames", {
    pred_kernel <- getPredKernel(params)
    expect_identical(dimnames(pred_kernel)$sp, 
                     dimnames(params@initial_n)$sp)
    expect_identical(dimnames(pred_kernel)$w_pred, 
                     dimnames(params@initial_n)$w)
    expect_identical(dimnames(pred_kernel)$w_prey, 
                     as.character(signif(params@w_full, 3)))
    expect_identical(pred_kernel(params), pred_kernel)
})
test_that("getting and setting pred kernel leads to same dynamics", {
    params <- NS_params_small
    params <- setPredKernel(params, pred_kernel = getPredKernel(params))
    sim1 <- project(NS_params_small, t_max = 0.1)
    sim2 <- project(params, t_max = 0.1)
    expect_equal(finalN(sim1), finalN(sim2), tolerance = 1e-4, ignore_attr = "params")
})

test_that("Can get and set pred_kernel slot", {
    params <- NS_params_small
    new <- 2 * pred_kernel(params)
    comment(new) <- "test"
    pred_kernel(params) <- new
    expect_identical(pred_kernel(params), new)
    expect_identical(getPredKernel(params), new)
})

## get_phi ----
test_that("get_phi works", {
    sp <- NS_species_params_small
    sp$pred_kernel_type <- "box"
    sp$ppmr_min <- 2
    sp$ppmr_max <- 4
    phi <- get_phi(sp, 1:5)
    expect_identical(phi[1, ], phi[2, ])
    expect_identical(phi[1, 1], 0)
    expect_identical(phi[1, 2], 1)
    expect_identical(phi[1, 5], 0)
    # call with invalid parameters
    sp$ppmr_max <- 1
    expect_error(get_phi(sp, 1:5),
                 "ppmr_min not less than ppmr_max")
    sp$pred_kernel_type <- "lognormal"
    sp$sigma <- 0
    expect_error(get_phi(sp, 1:5),
                 "The function lognormal_pred_kernel returned a zero predation kernel")
})

test_that("get_phi throws error if predation kernel parameters are missing", {
    sp <- NS_species_params_small
    sp$pred_kernel_type <- "box"
    # parameters missing entirely
    expect_error(get_phi(sp, 1:5),
                 "missing from the parameter dataframe: ppmr_min")
    # some entries missing
    sp$ppmr_min <- 2
    sp$ppmr_max <- 4
    sp$ppmr_min[2] <- NA
    expect_error(get_phi(sp, 1:5),
                 "arguments for the predation kernel function box_pred_kernel are NA")
})

test_that("default_pred_kernel_params sets defaults for data frames and params", {
    sp <- data.frame(species = c("A", "B"),
                     w_max = c(10, 20),
                     stringsAsFactors = FALSE)
    sp_defaulted <- default_pred_kernel_params(sp)
    expect_identical(sp_defaulted$pred_kernel_type, c("lognormal", "lognormal"))
    expect_identical(sp_defaulted$beta, c(30, 30))
    expect_identical(sp_defaulted$sigma, c(2, 2))

    params <- NS_params_small
    params@species_params$pred_kernel_type <- NA
    params@species_params$beta <- NA
    params@species_params$sigma <- NA
    params_defaulted <- default_pred_kernel_params(params)
    expect_true(all(params_defaulted@species_params$pred_kernel_type == "lognormal"))
    expect_true(all(params_defaulted@species_params$beta == 30))
    expect_true(all(params_defaulted@species_params$sigma == 2))
})

test_that("default_pred_kernel_params leaves manually set full kernels unchanged", {
    params <- setPredKernel(NS_params_small, pred_kernel = getPredKernel(NS_params_small))
    params@species_params$beta[] <- NA

    unchanged <- default_pred_kernel_params(params)

    expect_true(all(is.na(unchanged@species_params$beta)))
    expect_identical(unchanged@pred_kernel, params@pred_kernel)
})

## Higher-order (bin-integrated) quadrature ----
test_that("default scheme is first-order and leaves the slot untouched", {
    params <- NS_params_small
    # The default (bin_average = FALSE) leaves the Fourier kernels and the
    # second_order_w slot untouched
    expect_unchanged(setPredKernel(params), params)
    expect_false(setPredKernel(params)@second_order_w[["bin_average"]])
})

test_that("bin_average = TRUE changes the kernels", {
    params <- NS_params_small
    # setPredKernel() reads the bin_average flag from the slot
    params@second_order_w[["bin_average"]] <- TRUE
    p_hi <- setPredKernel(params)
    expect_false(isTRUE(all.equal(p_hi@ft_pred_kernel_e,
                                  NS_params_small@ft_pred_kernel_e)))
    expect_false(isTRUE(all.equal(p_hi@ft_pred_kernel_p,
                                  NS_params_small@ft_pred_kernel_p)))
    # The rates remain finite and a projection runs
    expect_true(all(is.finite(getEncounter(p_hi))))
    expect_true(all(is.finite(getPredRate(p_hi))))
    sim <- project(p_hi, t_max = 1, t_save = 1)
    expect_true(all(is.finite(sim@n)))
})

test_that("the high-order kernels persist through recalculation", {
    params <- NS_params_small
    # Enabling via the slot setter re-runs setParams() and builds the
    # high-order kernels
    second_order_w(params) <- c(bin_average = TRUE)
    hi_e <- params@ft_pred_kernel_e
    # A bare setPredKernel() recalculation keeps the high-order kernels
    expect_equal(setPredKernel(params)@ft_pred_kernel_e, hi_e)
    # Changing a kernel parameter re-runs setPredKernel() via the setParams
    # pipeline; the high-order scheme must be retained
    p2 <- params
    species_params(p2)$beta <- species_params(p2)$beta * 1.1
    expect_true(p2@second_order_w[["bin_average"]])
    p2_lo <- NS_params_small
    species_params(p2_lo)$beta <- species_params(p2_lo)$beta * 1.1
    expect_false(isTRUE(all.equal(p2@ft_pred_kernel_e, p2_lo@ft_pred_kernel_e)))
    # Switching back to first order also persists
    second_order_w(p2) <- c(bin_average = FALSE)
    expect_equal(p2@ft_pred_kernel_e, p2_lo@ft_pred_kernel_e)
})

test_that("bin-integrated box kernel matches the analytic weights", {
    params <- NS_params_small
    # Set the box parameters via the slot directly so that the kernel is not
    # recalculated before ppmr_min/ppmr_max are in place.
    params@species_params$pred_kernel_type <- "box"
    params@species_params$ppmr_min <- 1
    params@species_params$ppmr_max <- 1e8  # wide support: whole bins lie inside
    params@second_order_w[["bin_average"]] <- TRUE
    p_hi <- setPredKernel(params)
    no_w_full <- length(params@w_full)
    beta_grid <- params@w_full[2] / params@w_full[1]
    # Invert the FFT of the first species to recover the kernel weights
    phi_e <- Re(fft(p_hi@ft_pred_kernel_e[1, ], inverse = TRUE)) / no_w_full
    phi_p <- Re(fft(p_hi@ft_pred_kernel_p[1, ], inverse = TRUE)) / no_w_full
    phi_d <- Re(fft(p_hi@ft_pred_kernel_d[1, ], inverse = TRUE)) / no_w_full
    # Offset m = 0 (own bin) carries no encounter; a mid-range offset lies fully
    # inside the box support
    expect_equal(unname(phi_e[1]), 0, tolerance = 1e-8)
    expect_equal(unname(phi_e[5]), (beta_grid + 1) / 2, tolerance = 1e-4)
    # The predation weight for a fully-covered bin is exactly 1. phi_p is the
    # reversed kernel, so offset m = 1 sits in the last entry.
    expect_equal(unname(phi_p[no_w_full]), 1, tolerance = 1e-4)
    # The diffusion integrand carries one more power of prey size than the
    # encounter (w_p^2 dw_p), so a fully-covered bin gets the bin average of
    # (w_p/w_k)^2, which is (beta^3 - 1)/(3(beta - 1)) = (beta^2 + beta + 1)/3,
    # the diffusion analogue of the encounter's (beta + 1)/2.
    expect_equal(unname(phi_d[5]),
                 (beta_grid^2 + beta_grid + 1) / 3, tolerance = 1e-4)
})

test_that("ft_pred_kernel_d equals ft_pred_kernel_e in the first-order scheme", {
    # With bin_average off the diffusion kernel mirrors the encounter kernel, so
    # the predation-diffusion integral is byte-identical to previous mizer.
    p <- setPredKernel(NS_params_small)
    expect_identical(p@ft_pred_kernel_d, p@ft_pred_kernel_e)
    # And it differs once bin-averaging is on (the extra power of prey size).
    p_hi <- NS_params_small
    p_hi@second_order_w[["bin_average"]] <- TRUE
    p_hi <- setPredKernel(p_hi)
    expect_false(isTRUE(all.equal(p_hi@ft_pred_kernel_d, p_hi@ft_pred_kernel_e)))
})
