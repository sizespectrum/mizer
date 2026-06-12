# second_order_w ----
test_that("second_order_w slot defaults to c(flux_limiter=FALSE, bin_average=FALSE)", {
    params <- NS_params_small
    sow <- second_order_w(params)
    expect_identical(sow, c(flux_limiter = FALSE, bin_average = FALSE))
})

test_that("second_order_w getter works", {
    params <- NS_params_small
    params@second_order_w <- c(flux_limiter = TRUE, bin_average = TRUE)
    sow <- second_order_w(params)
    expect_identical(sow, c(flux_limiter = TRUE, bin_average = TRUE))
})

test_that("second_order_w setter with single TRUE sets both entries", {
    params <- NS_params_small
    second_order_w(params) <- TRUE
    expect_identical(second_order_w(params),
                     c(flux_limiter = TRUE, bin_average = TRUE))
})

test_that("second_order_w setter with single FALSE sets both entries", {
    params <- NS_params_small
    second_order_w(params) <- TRUE
    second_order_w(params) <- FALSE
    expect_identical(second_order_w(params),
                     c(flux_limiter = FALSE, bin_average = FALSE))
})

test_that("second_order_w setter with named vector sets individual entries", {
    params <- NS_params_small
    second_order_w(params) <- c(flux_limiter = TRUE)
    expect_true(second_order_w(params)[["flux_limiter"]])
    expect_false(second_order_w(params)[["bin_average"]])

    second_order_w(params) <- c(bin_average = TRUE)
    expect_true(second_order_w(params)[["flux_limiter"]])
    expect_true(second_order_w(params)[["bin_average"]])
})

test_that("second_order_w setter validates input", {
    params <- NS_params_small
    expect_error(second_order_w(params) <- "yes")
    expect_error(second_order_w(params) <- NA)
    expect_error(second_order_w(params) <- c(unknown = TRUE))
})

test_that("second_order_w setter re-runs setParams", {
    params <- NS_params_small
    second_order_w(params) <- TRUE
    expect_identical(second_order_w(params),
                     c(flux_limiter = TRUE, bin_average = TRUE))
    expect_s4_class(params, "MizerParams")
})

test_that("flux_limiter can be set to a reconstruction scheme name", {
    params <- NS_params_small
    # TRUE / van_leer / centred all switch the flag on; only the name persisted
    # in metadata distinguishes the reconstruction.
    second_order_w(params) <- c(flux_limiter = TRUE)
    expect_identical(mizer:::flux_limiter_scheme(params), "van_leer")

    second_order_w(params) <- c(flux_limiter = "centred")
    expect_true(second_order_w(params)[["flux_limiter"]])
    expect_identical(mizer:::flux_limiter_scheme(params), "centred")
    expect_identical(params@metadata$flux_reconstruction, "centred")

    # A bare scheme name sets only the flux_limiter entry.
    p2 <- NS_params_small
    second_order_w(p2) <- "centred"
    expect_identical(mizer:::flux_limiter_scheme(p2), "centred")
    expect_false(second_order_w(p2)[["bin_average"]])

    # Switching off clears the stored reconstruction.
    second_order_w(params) <- c(flux_limiter = "none")
    expect_false(second_order_w(params)[["flux_limiter"]])
    expect_identical(mizer:::flux_limiter_scheme(params), "none")
    expect_null(params@metadata$flux_reconstruction)

    expect_error(second_order_w(params) <- c(flux_limiter = "bogus"))
})

test_that("flux_reconstruction persists when bin_average rebuilds the model", {
    params <- NS_params_small
    second_order_w(params) <- c(flux_limiter = "centred", bin_average = TRUE)
    expect_identical(mizer:::flux_limiter_scheme(params), "centred")
    expect_true(second_order_w(params)[["bin_average"]])
})

test_that("setting bin_average rebuilds the higher-order predation kernels", {
    params <- NS_params_small
    p_hi <- params
    second_order_w(p_hi) <- c(bin_average = TRUE)
    # The bin-averaged quadrature gives different Fourier kernels
    expect_false(isTRUE(all.equal(p_hi@ft_pred_kernel_e,
                                  params@ft_pred_kernel_e)))
    # It agrees with building the kernel from the slot directly
    p_direct <- params
    p_direct@second_order_w[["bin_average"]] <- TRUE
    expect_equal(p_hi@ft_pred_kernel_e,
                 setPredKernel(p_direct)@ft_pred_kernel_e)
    # Toggling back restores the first-order kernels
    second_order_w(p_hi) <- c(bin_average = FALSE)
    expect_equal(p_hi@ft_pred_kernel_e, params@ft_pred_kernel_e)
})

test_that("second_order_w slot is preserved by upgradeParams", {
    params <- NS_params_small
    params2 <- upgradeParams(params)
    expect_identical(second_order_w(params2),
                     c(flux_limiter = FALSE, bin_average = FALSE))
})

test_that("validObject accepts second_order_w = c(TRUE, TRUE)", {
    params <- NS_params_small
    params@second_order_w <- c(flux_limiter = TRUE, bin_average = TRUE)
    expect_silent(validObject(params))
})

test_that("validObject rejects invalid second_order_w", {
    params <- NS_params_small
    params@second_order_w <- logical(0)
    expect_error(validObject(params))

    params@second_order_w <- TRUE  # unnamed single value
    expect_error(validObject(params))
})
