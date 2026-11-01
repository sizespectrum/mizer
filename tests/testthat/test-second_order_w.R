# second_order_w ----
test_that("second_order_w slot defaults to c(flux='upwind', bin_average='FALSE')", {
    params <- NS_params_small
    sow <- second_order_w(params)
    expect_identical(sow, list(flux = "upwind", bin_average = FALSE))
})

test_that("second_order_w getter works", {
    params <- NS_params_small
    params@second_order_w <- list(flux = "van_leer", bin_average = TRUE)
    sow <- second_order_w(params)
    expect_identical(sow, list(flux = "van_leer", bin_average = TRUE))
})

test_that("second_order_w setter with single TRUE sets both entries", {
    params <- NS_params_small
    second_order_w(params) <- TRUE
    expect_identical(second_order_w(params),
                     list(flux = "van_leer", bin_average = TRUE))
})

test_that("second_order_w setter with single FALSE sets both entries", {
    params <- NS_params_small
    second_order_w(params) <- TRUE
    second_order_w(params) <- FALSE
    expect_identical(second_order_w(params),
                     list(flux = "upwind", bin_average = FALSE))
})

test_that("second_order_w setter with named vector sets individual entries", {
    params <- NS_params_small
    second_order_w(params) <- c(flux = "van_leer")
    expect_identical(second_order_w(params)[["flux"]], "van_leer")
    expect_false(second_order_w(params)[["bin_average"]])

    second_order_w(params) <- c(bin_average = TRUE)
    expect_identical(second_order_w(params)[["flux"]], "van_leer")
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
                     list(flux = "van_leer", bin_average = TRUE))
    expect_s3_class(params, "MizerParams")
})

test_that("flux can be set to a reconstruction scheme name", {
    params <- NS_params_small
    second_order_w(params) <- c(flux = TRUE)
    expect_identical(mizer:::flux_limiter_scheme(params), "van_leer")

    second_order_w(params) <- c(flux = "centred")
    expect_identical(second_order_w(params)[["flux"]], "centred")
    expect_identical(mizer:::flux_limiter_scheme(params), "centred")
    expect_null(params@metadata$flux_reconstruction)

    # A bare scheme name sets only the flux entry.
    p2 <- NS_params_small
    second_order_w(p2) <- "centred"
    expect_identical(mizer:::flux_limiter_scheme(p2), "centred")
    expect_false(second_order_w(p2)[["bin_average"]])

    # Switching off reverts to upwind.
    second_order_w(params) <- c(flux = "upwind")
    expect_identical(second_order_w(params)[["flux"]], "upwind")
    expect_identical(mizer:::flux_limiter_scheme(params), "none")

    expect_error(second_order_w(params) <- c(flux = "bogus"))
})

test_that("flux scheme persists when bin_average rebuilds the model", {
    params <- NS_params_small
    second_order_w(params) <- c(flux = "centred", bin_average = TRUE)
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

test_that("default_second_order_w gives the scheme of the defaults edition", {
    first <- list(flux = "upwind", bin_average = FALSE)
    second <- list(flux = "van_leer", bin_average = TRUE)
    expect_identical(default_second_order_w(1), first)
    expect_identical(default_second_order_w(2), first)
    expect_identical(default_second_order_w(3), second)
    expect_identical(default_second_order_w(4), second)
})

test_that("resolve_second_order_w resolves against the edition's default", {
    # NULL gives the default scheme of the edition
    expect_identical(resolve_second_order_w(NULL, 2),
                     list(flux = "upwind", bin_average = FALSE))
    expect_identical(resolve_second_order_w(NULL, 3),
                     list(flux = "van_leer", bin_average = TRUE))
    # A single logical sets both entries in any edition
    expect_identical(resolve_second_order_w(FALSE, 3),
                     list(flux = "upwind", bin_average = FALSE))
    expect_identical(resolve_second_order_w(TRUE, 1),
                     list(flux = "van_leer", bin_average = TRUE))
    # A value that sets one entry leaves the other at the edition's default
    expect_identical(resolve_second_order_w("centred", 3),
                     list(flux = "centred", bin_average = TRUE))
    expect_identical(resolve_second_order_w("centred", 1),
                     list(flux = "centred", bin_average = FALSE))
    expect_identical(resolve_second_order_w(c(bin_average = FALSE), 3),
                     list(flux = "van_leer", bin_average = FALSE))
    expect_error(resolve_second_order_w("bogus", 3))
})

test_that("upgradeParams is idempotent for current-format second_order_w", {
    params <- NS_params_small
    params2 <- upgradeParams(params)
    expect_identical(second_order_w(params), second_order_w(params2))
})

test_that("validObject accepts valid second_order_w", {
    params <- NS_params_small
    params@second_order_w <- list(flux = "van_leer", bin_average = TRUE)
    expect_silent(validObject(params))
})

test_that("validObject rejects invalid second_order_w", {
    params <- NS_params_small
    params@second_order_w <- list()
    expect_error(validObject(params))

    params@second_order_w <- list(flux = "upwind", bin_average = NA)
    expect_error(validObject(params))

    params@second_order_w <- list(flux = "bogus", bin_average = FALSE)
    expect_error(validObject(params))
})
