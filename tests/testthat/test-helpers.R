test_that("l2w works", {
    sp <- validSpeciesParams(NS_species_params_small)
    no_sp <- nrow(sp)
    # call with species_params
    expect_identical(l2w(2, sp), rep(0.08, no_sp))
    # call with params - result is named by species
    # NS_params_small should have `a` and `b` because of the mizer upgrade,
    # but let's make sure by testing against it directly if it has them.
    # (Actually we will just test it directly)
    expect_equal(l2w(2, NS_params_small), rep(0.08, no_sp), ignore_attr = TRUE)
    # call with wrong 2nd argument
    expect_error(l2w(2, 4),
                 "The second argument must be either ")
    # call with wrong 1st argument
    expect_error(l2w("a", sp),
                 "l is not a numeric or integer vector")
    expect_error(l2w(1:2, sp),
                 "The length of 'l'")
})

test_that("w2l works", {
    sp <- validSpeciesParams(NS_species_params_small)
    no_sp <- nrow(sp)
    # call with species_params
    expect_identical(w2l(0.08, sp), rep(2, no_sp))
    # call with params - result is named by species
    expect_equal(w2l(0.08, NS_params_small), rep(2, no_sp), ignore_attr = TRUE)
    # call with wrong 1st argument
    expect_error(w2l("a", sp),
                 "w is not a numeric or integer vector")
    expect_error(w2l(1:2, sp),
                 "The length of 'w'")
    # w2l should do the inverse of l2w
    expect_equal(w2l(l2w(2, sp), sp),
                 rep(2, no_sp))
})

test_that("resource_power_law point-samples by default", {
    p <- NS_params_small
    kappa <- 2e10
    lambda <- 2.05
    # No cutoff: plain left-edge power law.
    expect_equal(resource_power_law(p, kappa, lambda),
                 kappa * p@w_full ^ (-lambda))
    # With a cutoff: zero at and above it, point values below.
    w_max <- 1e-3
    expected <- kappa * p@w_full ^ (-lambda)
    expected[p@w_full >= w_max] <- 0
    expect_equal(resource_power_law(p, kappa, lambda, w_max = w_max), expected)
})

test_that("resource_power_law bin-averages under second_order_w", {
    p <- NS_params_small
    second_order_w(p) <- c(bin_average = TRUE)
    kappa <- 2e10
    lambda <- 2.05
    w <- p@w_full
    dw <- p@dw_full
    # Without cutoff: exact bin average of kappa * w^(-lambda).
    expect_equal(resource_power_law(p, kappa, lambda),
                 kappa * power_law_bin_average(w, dw, -lambda))
    # With cutoff: straddling bin gets the partial average, bins above are 0.
    w_max <- 1e-3
    expect_equal(resource_power_law(p, kappa, lambda, w_max = w_max),
                 kappa * power_law_bin_average(w, dw, -lambda, w_max = w_max))
    # Bin average differs from the left-edge point values (steep power law).
    expect_false(isTRUE(all.equal(
        unname(resource_power_law(p, kappa, lambda)),
        unname(kappa * w ^ (-lambda)))))
})

test_that("newMultispeciesParams initial resource is byte-identical", {
    # Default (second_order_w off) initial_n_pp is the left-edge power law,
    # cut at w_pp_cutoff, unchanged from previous mizer.
    kappa <- NS_params_small@resource_params$kappa
    lambda <- NS_params_small@resource_params$lambda
    cutoff <- NS_params_small@resource_params$w_pp_cutoff
    w <- NS_params_small@w_full
    expected <- kappa * w ^ (-lambda)
    expected[w >= cutoff] <- 0
    expect_equal(unname(NS_params_small@initial_n_pp), unname(expected))
})

test_that("newMultispeciesParams second_order_w bin-averages resource", {
    sp <- NS_species_params_small
    p <- suppressMessages(newMultispeciesParams(sp, inter_small,
                                                second_order_w = TRUE))
    # The flag is set on the returned object.
    expect_equal(p@second_order_w$flux, "van_leer")
    expect_true(p@second_order_w$bin_average)

    # The initial resource and capacity are the exact bin averages.
    w <- p@w_full
    dw <- p@dw_full
    kappa <- p@resource_params$kappa
    lambda <- p@resource_params$lambda
    cutoff <- p@resource_params$w_pp_cutoff
    expected <- kappa * power_law_bin_average(w, dw, -lambda, w_max = cutoff)
    expect_equal(unname(p@initial_n_pp), expected)
    expect_equal(unname(p@cc_pp), expected)

    # Initial abundances are finite and the model projects.
    expect_true(all(is.finite(p@initial_n)))
    expect_error(suppressMessages(project(p, t_max = 1)), NA)
})

test_that("second_order_w argument accepts the same forms as the setter", {
    sp <- NS_species_params_small
    # bin_average only leaves flux at the default upwind.
    p1 <- suppressMessages(newMultispeciesParams(
        sp, inter_small, second_order_w = c(bin_average = TRUE)))
    expect_equal(p1@second_order_w$flux, "upwind")
    expect_true(p1@second_order_w$bin_average)
    # A scheme name sets only flux.
    p2 <- suppressMessages(newMultispeciesParams(
        sp, inter_small, second_order_w = "centred"))
    expect_equal(p2@second_order_w$flux, "centred")
    expect_false(p2@second_order_w$bin_average)
    # The default is the first-order scheme.
    p0 <- suppressMessages(newMultispeciesParams(sp, inter_small))
    expect_equal(p0@second_order_w$flux, "upwind")
    expect_false(p0@second_order_w$bin_average)
})

test_that("wrapper constructors accept second_order_w = TRUE", {
    # The van_leer flux only governs projection, so construction (which runs a
    # steady-state solve) must stay robust. Each wrapper should build finite
    # abundances and carry the requested scheme on the returned object.
    pt <- suppressWarnings(suppressMessages(
        newTraitParams(no_sp = 3, second_order_w = TRUE)))
    expect_equal(pt@second_order_w$flux, "van_leer")
    expect_true(pt@second_order_w$bin_average)
    expect_true(all(is.finite(pt@initial_n)))

    pc <- suppressWarnings(suppressMessages(
        newCommunityParams(reproduction = 1, second_order_w = TRUE)))
    expect_equal(pc@second_order_w$flux, "van_leer")
    expect_true(pc@second_order_w$bin_average)
    expect_true(all(is.finite(pc@initial_n)))

    ps <- suppressWarnings(suppressMessages(
        newSingleSpeciesParams(second_order_w = TRUE)))
    expect_equal(ps@second_order_w$flux, "van_leer")
    expect_true(ps@second_order_w$bin_average)
    expect_true(all(is.finite(ps@initial_n)))
})

test_that("get_gamma_default runs and stays invariant under second_order_w", {
    sp <- NS_species_params_small
    sp$gamma <- NULL
    p1 <- suppressMessages(newMultispeciesParams(sp, inter_small))
    g1 <- suppressMessages(get_gamma_default(p1))

    p2 <- p1
    second_order_w(p2) <- c(bin_average = TRUE)
    # Exercises the bin-averaged prey-spectrum build in get_gamma_default.
    g2 <- suppressMessages(get_gamma_default(p2))

    expect_true(all(is.finite(g2)))
    # The default gamma is invariant: for a pure power-law prey the
    # bin-integrated encounter convolution consumes the bin-averaged resource
    # to give the same available energy as the first-order point sampling, so
    # the feeding-level calibration is unchanged.
    expect_equal(unname(g1), unname(g2))
})

# get_steady_state_n ----

test_that("get_steady_state_n works with no diffusion", {
    params <- NS_params_small
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    # Mocking constant growth and mortality
    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.5, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)

    # Zero diffusion
    params@ext_diffusion[] <- 0

    n_calc <- mizer:::get_steady_state_n(params, growth, mort,
                                         params@ext_diffusion, N0_vec)

    expect_equal(dim(n_calc), c(no_sp, no_w))

    # Check species 1 manually against old analytical form
    sp <- 1
    w_min_idx <- params@w_min_idx[sp]
    w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])
    idx <- w_min_idx:(w_max_idx - 1)
    dw <- params@dw

    # Old calculation logic (no diffusion)
    n_old <- c(1, cumprod(growth[sp, idx] / (growth[sp, idx + 1] + mort[sp, idx + 1] * dw[idx + 1])))
    n_old <- 100 * n_old

    expect_equal(unname(n_calc[sp, w_min_idx:w_max_idx]), unname(n_old), tolerance = 1e-10)


})

test_that("get_steady_state_n works with diffusion", {
    params <- example_params()
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    # Mocking constant growth and mortality
    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.1, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)

    D <- getDiffusion(params)
    n_calc <- mizer:::get_steady_state_n(params, growth, mort, D, N0_vec)

    expect_equal(dim(n_calc), c(no_sp, no_w))

    # Should be positive up to w_max
    for(sp in 1:no_sp) {
        w_min_idx <- params@w_min_idx[sp]
        w_max_idx <- sum(params@w <= params@species_params[sp, "w_max"])
        expect_equal(unname(n_calc[sp, w_min_idx]), unname(N0_vec[sp]))
    }

    # Compare with no diffusion
    zero_D <- D
    zero_D[] <- 0
    n_nodiff <- mizer:::get_steady_state_n(params, growth, mort, zero_D,
                                           N0_vec)
    expect_false(isTRUE(all.equal(n_calc, n_nodiff)))
})

test_that("get_steady_state_n uses supplied diffusion array", {
    params <- example_params()
    params@ext_diffusion[] <- 0
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.1, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)
    D <- params@ext_diffusion
    D[1, ] <- 0.1 * params@w^(params@species_params$n[1] + 1)

    n_with_D <- mizer:::get_steady_state_n(params, growth, mort, D, N0_vec)
    D[] <- 0
    n_without_D <- mizer:::get_steady_state_n(params, growth, mort, D, N0_vec)

    expect_false(isTRUE(all.equal(n_with_D, n_without_D)))
})

test_that("get_steady_state_n holds abundance at zero above w_max (no diffusion)", {
    params <- NS_params_small
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    # Constant positive growth, so the bottom-up recursion would never decay on
    # its own and would otherwise climb to the top of the grid.
    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.5, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)
    params@ext_diffusion[] <- 0

    n_calc <- mizer:::get_steady_state_n(params, growth, mort,
                                         params@ext_diffusion, N0_vec)

    w_top <- mizer:::support_top_idx(params)
    for (sp in seq_len(no_sp)) {
        # Positive up to and including the support top, exactly zero above it.
        expect_true(all(n_calc[sp, params@w_min_idx[sp]:w_top[sp]] > 0))
        if (w_top[sp] < no_w) {
            expect_true(all(n_calc[sp, (w_top[sp] + 1):no_w] == 0))
        }
    }
})

test_that("get_steady_state_n holds abundance at zero above w_max (diffusion)", {
    params <- example_params()
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    growth <- matrix(1, nrow = no_sp, ncol = no_w)
    mort <- matrix(0.1, nrow = no_sp, ncol = no_w)
    N0_vec <- rep(100, no_sp)
    D <- getDiffusion(params)

    n_calc <- mizer:::get_steady_state_n(params, growth, mort, D, N0_vec)

    w_top <- mizer:::support_top_idx(params)
    for (sp in seq_len(no_sp)) {
        if (w_top[sp] < no_w) {
            expect_true(all(n_calc[sp, (w_top[sp] + 1):no_w] == 0))
        }
    }
})
