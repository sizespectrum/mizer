# getDiffusion ----
test_that("getDiffusion returns ArraySpeciesBySize with correct dimensions", {
    params <- single_sp_params
    d <- getDiffusion(params)
    expect_true(is.ArraySpeciesBySize(d))
    expect_identical(dim(d), dim(params@initial_n))
})

test_that("getDiffusion includes ext_diffusion", {
    params <- single_sp_params
    d_base <- getDiffusion(params)
    # Adding a constant to ext_diffusion should shift getDiffusion by the same amount
    params@ext_diffusion[] <- 1
    d_with_ext <- getDiffusion(params)
    expect_equal(d_with_ext, d_base + 1, ignore_attr = TRUE)
})

test_that("getDiffusion dispatches via rates_funcs", {
    params <- single_sp_params
    e <- globalenv()
    e$constant_diffusion <- function(params, n, n_pp, n_other, t, feeding_level, ...) {
        array(42, dim = dim(params@initial_n), dimnames = dimnames(params@initial_n))
    }
    params@rates_funcs$Diffusion <- "constant_diffusion"
    d <- getDiffusion(params)
    expect_true(all(d == 42))
})

test_that("r$diffusion is included in getRates output", {
    params <- single_sp_params
    r <- getRates(params)
    expect_true("diffusion" %in% names(r))
    expect_identical(dim(r$diffusion), dim(params@initial_n))
})

test_that("r$diffusion matches getDiffusion", {
    params <- single_sp_params
    r <- getRates(params)
    expect_equal(r$diffusion, getDiffusion(params), ignore_attr = TRUE)
})

# mizerDiffusion behaviour ----

test_that("mizerDiffusion accepts pre-computed feeding_level and gives same result", {
    params <- single_sp_params
    params@use_predation_diffusion <- TRUE
    n    <- initialN(params)
    n_pp <- initialNResource(params)
    fl   <- getFeedingLevel(params, n = n, n_pp = n_pp)
    d_auto  <- getDiffusion(params, n, n_pp)
    d_given <- mizerDiffusion(params, n = n, n_pp = n_pp,
                               n_other = initialNOther(params),
                               t = 0, feeding_level = fl)
    expect_equal(d_auto, d_given, ignore_attr = TRUE)
})

test_that("mizerDiffusion is zero when feeding level is 1 everywhere", {
    # D(w) = (1 - f(w)) * ... so f = 1 gives D = ext_diffusion
    params <- single_sp_params
    params@use_predation_diffusion <- TRUE
    n    <- initialN(params)
    n_pp <- initialNResource(params)
    fl_one <- matrix(1, nrow = nrow(n), ncol = ncol(n), dimnames = dimnames(n))
    d <- mizerDiffusion(params, n = n, n_pp = n_pp,
                         n_other = initialNOther(params),
                         t = 0, feeding_level = fl_one)
    expect_equal(d, params@ext_diffusion, ignore_attr = TRUE)
})

test_that("mizerDiffusion scales as alpha^2", {
    params <- single_sp_params
    params@use_predation_diffusion <- TRUE
    n    <- initialN(params)
    n_pp <- initialNResource(params)
    d1 <- getDiffusion(params, n, n_pp)
    # Double alpha → diffusion should quadruple (alpha enters as alpha^2)
    params2 <- params
    params2@species_params$alpha <- params@species_params$alpha * 2
    # Keep search_vol unchanged to isolate the alpha^2 factor
    params2@search_vol[] <- params@search_vol
    d2 <- getDiffusion(params2, n, n_pp)
    expect_equal(d2, d1 * 4, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("use_predation_diffusion<- updates `time_modified`", {
    params <- single_sp_params
    p2 <- params
    use_predation_diffusion(p2) <- TRUE
    expect_false(identical(p2@time_modified, params@time_modified))
})

test_that("mizerDiffusion increases when fish prey are present", {
    # With use_predation_diffusion = TRUE, adding fish as prey increases D.
    # This exercises the params@interaction %*% n term.
    params <- newMultispeciesParams(NS_species_params_gears_small, inter_small,
                                    no_w = 30, info_level = 0)
    params@use_predation_diffusion <- TRUE
    n_zero <- initialN(params)
    n_zero[] <- 0
    n_pp <- initialNResource(params)
    d_no_fish   <- getDiffusion(params, n_zero, n_pp)
    d_with_fish <- getDiffusion(params, initialN(params), n_pp)
    # At least one species should have larger diffusion with fish present
    expect_true(any(d_with_fish > d_no_fish))
})

test_that("predation diffusion uses its own bin-integrated kernel (#384)", {
    params <- single_sp_params
    params@use_predation_diffusion <- TRUE
    n <- initialN(params)
    n_pp <- initialNResource(params)
    d_lo <- getDiffusion(params, n, n_pp)

    p_hi <- params
    second_order_w(p_hi) <- c(bin_average = TRUE)
    d_hi <- getDiffusion(p_hi, n, n_pp)

    # The dedicated diffusion kernel (beta^{3s} Jacobian) shifts the
    # predation-diffusion rate away from the first-order point-sampled value.
    expect_true(all(is.finite(d_hi)))
    expect_true(all(d_hi >= 0))
    expect_false(isTRUE(all.equal(unname(d_hi), unname(d_lo))))
})

test_that("getDiffusion works with a custom predation kernel (#373)", {
    # When the user has set a custom pred_kernel (with a comment) that does not
    # depend only on the predator/prey size ratio, the FFT method cannot be
    # used. getDiffusion should then use the same direct summation over the
    # full predation kernel as projectEncounter. We check the result against a
    # brute-force evaluation of the diffusion integral.
    params <- newMultispeciesParams(NS_species_params_gears_small, inter_small,
                                    no_w = 30, info_level = 0)
    params@use_predation_diffusion <- TRUE
    # Make the kernel genuinely depend on predator size (not just the ratio) by
    # scaling each predator species' kernel, then store it as a custom kernel.
    pk <- getPredKernel(params)
    pk["Herring", , ] <- pk["Herring", , ] * 1.5
    params <- setPredKernel(params, pred_kernel = pk)
    expect_false(is.null(comment(params@pred_kernel)))

    n <- initialN(params)
    n_pp <- initialNResource(params)
    d <- getDiffusion(params, n, n_pp)

    # Brute-force diffusion integral I_d(w) = int phi(w, w_p) N(w_p) w_p^2 dw_p
    idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
    no_sp <- nrow(n)
    integral_d <- matrix(0, nrow = no_sp, ncol = length(params@w))
    for (i in seq_len(no_sp)) {
        prey <- params@species_params$interaction_resource[i] * n_pp *
            params@w_full^2 * params@dw_full
        prey[idx_sp] <- prey[idx_sp] +
            (params@interaction[i, ] %*% n) * params@w^2 * params@dw
        for (k in seq_along(params@w)) {
            integral_d[i, k] <- sum(params@pred_kernel[i, k, ] * prey)
        }
    }
    fl <- getFeedingLevel(params, n = n, n_pp = n_pp)
    expected <- (1 - fl) * params@search_vol *
        params@species_params$alpha^2 * integral_d + params@ext_diffusion

    expect_equal(d, expected, ignore_attr = TRUE)
})

test_that("predation diffusion default path mirrors the encounter kernel (#384)", {
    # With bin_average off, ft_pred_kernel_d == ft_pred_kernel_e, so the
    # diffusion integral is unchanged from previous mizer.
    params <- single_sp_params
    params@use_predation_diffusion <- TRUE
    expect_identical(params@ft_pred_kernel_d, params@ft_pred_kernel_e)
    n <- initialN(params)
    n_pp <- initialNResource(params)
    expect_true(all(is.finite(getDiffusion(params, n, n_pp))))
})
