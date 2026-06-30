# Tests for the direct (Newton) steady-state solver.

# A steadied model is a good starting guess for the root finder.
p_steady <- steady(NS_params_small, t_max = 50, progress_bar = FALSE)

test_that("steadyNewton returns a fixed point of the dynamics", {
    pn <- steadyNewton(p_steady)
    expect_s4_class(pn, "MizerParams")

    # The returned spectra should barely move when projected forward.
    sim <- project(pn, t_max = 1, dt = 0.25, t_save = 1)
    n0 <- pn@initial_n
    n1 <- finalN(sim)
    support <- n0 > 0
    drift <- max((abs(n1 - n0) / n0)[support])
    expect_lt(drift, 1e-3)
})

test_that("steadyNewton solves the discrete steady-state equation to tolerance", {
    # Solve the held-RDD problem directly and check the residual at the root.
    params <- validParams(p_steady)
    rdd_const <- getRDD(params)
    active <- mizer:::steady_active_set(params)
    resfn <- mizer:::steady_state_residual(params, rdd_const,
                                           params@initial_n_other,
                                           params@initial_effort, active)
    x0 <- log(params@initial_n[active$mask])
    sol <- nleqslv::nleqslv(x0, resfn, method = "Newton", global = "dbldog",
                            control = list(maxit = 100, ftol = 1e-10,
                                           xtol = 1e-10))
    expect_lte(sol$termcd, 2)
    expect_lt(max(abs(sol$fvec)), 1e-8)
})

test_that("steadyNewton agrees with steady on a stable model", {
    pn <- steadyNewton(p_steady)
    # Compare where the steady() spectra are well above zero.
    # steadyNewton finds a tighter fixed point than steady() reaches by
    # time-stepping, so a few percent difference on the well-resolved bins is
    # expected (and is the point of the exercise).
    big <- p_steady@initial_n > max(p_steady@initial_n) * 1e-6
    rel <- abs(pn@initial_n - p_steady@initial_n) / p_steady@initial_n
    expect_lt(max(rel[big]), 0.1)
})

test_that("steadyNewton honours the preserve and reproduction arguments", {
    pn_level <- steadyNewton(p_steady, preserve = "reproduction_level")
    expect_equal(getReproductionLevel(pn_level),
                 getReproductionLevel(p_steady), tolerance = 1e-5)

    pn_rmax <- steadyNewton(p_steady, preserve = "R_max")
    expect_equal(pn_rmax@species_params$R_max,
                 p_steady@species_params$R_max, tolerance = 1e-5)

    pn_erepro <- suppressWarnings(steadyNewton(p_steady, preserve = "erepro"))
    expect_equal(pn_erepro@species_params$erepro,
                 p_steady@species_params$erepro, tolerance = 1e-5)

    pn_none <- steadyNewton(p_steady, reproduction = "dynamic")
    expect_equal(pn_none@species_params$R_max,
                 p_steady@species_params$R_max, tolerance = 1e-5)
    expect_equal(pn_none@species_params$erepro,
                 p_steady@species_params$erepro, tolerance = 1e-5)

    # Verify preserve argument is ignored when reproduction = "dynamic"
    pn_ignored <- steadyNewton(p_steady, reproduction = "dynamic", preserve = "invalid_option")
    expect_equal(pn_ignored@species_params$R_max,
                 p_steady@species_params$R_max, tolerance = 1e-5)

    # Verify verbose = TRUE captures iteration report output
    out <- capture.output(steadyNewton(p_steady, verbose = TRUE))
    expect_true(any(grepl("Iteration report", out)))
})

test_that("steadyNewton handles extinctions under reproduction = 'dynamic' with relative floor", {
    # Make species 3 (Cod) unviable by setting its reproduction efficiency extremely low
    p_extinct <- p_steady
    p_extinct@species_params$erepro[3] <- 1e-12
    p_extinct@initial_n[3, ] <- p_steady@initial_n[3, ]

    # Verify that the solver issues the extinction warning and pegs to the floor
    expect_warning(pn_ext <- steadyNewton(p_extinct, reproduction = "dynamic", extinction_floor = 1e-6),
                   "went extinct and were pegged to their abundance floor")

    # Verify Cod abundance is pegged exactly to the floor (1e-6 of its initial abundance)
    lo <- p_extinct@w_min_idx[3]
    ratio <- pn_ext@initial_n[3, lo] / p_extinct@initial_n[3, lo]
    expect_equal(ratio, 1e-6, tolerance = 1e-2)
})

test_that("steadyNewton errors for unsupported resource dynamics", {
    p_log <- setResource(NS_params_small,
                         resource_dynamics = "resource_logistic")
    expect_error(steadyNewton(p_log), "semichemostat")
})

test_that("steadyNewton works with the second-order (van Leer) scheme", {
    p <- NS_params_small
    sow <- second_order_w(p)
    sow$flux <- "van_leer"
    sow$bin_average <- TRUE
    second_order_w(p) <- sow
    ps <- steady(p, t_max = 100, progress_bar = FALSE)

    pn <- steadyNewton(ps)
    expect_s4_class(pn, "MizerParams")
    # The returned state must still be a fixed point of the (second-order)
    # dynamics. The van Leer limiter is only Lipschitz, so the Newton residual
    # cannot be driven to machine precision, but the projected drift is the
    # honest test and stays small.
    sim <- project(pn, t_max = 1, dt = 0.25, t_save = 1)
    n0 <- pn@initial_n
    n1 <- finalN(sim)
    support <- n0 > 0
    drift <- max((abs(n1 - n0) / n0)[support])
    expect_lt(drift, 1e-3)
})

test_that("support_top_idx drops the pile-up bin for the second-order scheme", {
    # First-order upwind feeds a class from the growth of the class below, so the
    # support reaches one bin past w_max; the second-order scheme feeds a class
    # from the growth at its own lower face, so it stops at w_max.
    p1 <- NS_params_small
    p2 <- NS_params_small
    sow <- second_order_w(p2)
    sow$flux <- "van_leer"
    second_order_w(p2) <- sow

    w_max_idx <- sapply(seq_len(nrow(p1@species_params)), function(i) {
        sum(p1@w <= p1@species_params$w_max[i])
    })
    no_w <- length(p1@w)
    expect_equal(mizer:::support_top_idx(p1), pmin(w_max_idx + 1L, no_w),
                 ignore_attr = TRUE)
    expect_equal(mizer:::support_top_idx(p2), pmin(w_max_idx, no_w),
                 ignore_attr = TRUE)
})

test_that("steady_active_set always reaches the grid truncation limit", {
    # Under the new dynamic support design, the active set always reaches
    # the grid truncation limit (support_top_idx) to allow the solver to
    # automatically discover the non-zero region.
    p <- NS_params_small
    no_w <- length(p@w)
    grid_top <- mizer:::support_top_idx(p)

    active0 <- mizer:::steady_active_set(p)
    expect_equal(active0$w_top, unname(grid_top))

    # Even if we zero out the tail of the abundances, the mask still reaches grid_top
    cutoff <- unname(pmax(p@w_min_idx + 2L, grid_top - 3L))
    for (i in seq_len(nrow(p@species_params))) {
        if (cutoff[i] < no_w) {
            p@initial_n[i, (cutoff[i] + 1L):no_w] <- 0
        }
    }
    active <- mizer:::steady_active_set(p)
    expect_equal(active$w_top, unname(grid_top))
})

test_that("support_top_idx is the first class above w_max", {
    p <- NS_params_small
    g <- getEGrowth(p)
    w_top <- mizer:::support_top_idx(p)
    no_w <- length(p@w)
    for (i in seq_len(nrow(p@species_params))) {
        # The support top is one class above w_max (the pile-up bin), capped at
        # the top of the grid.
        w_max_idx <- sum(p@w <= p@species_params$w_max[i])
        expect_equal(unname(w_top[i]), min(w_max_idx + 1, no_w))
        # For a standard model this coincides with the first zero-growth class.
        if (w_top[i] < no_w) {
            expect_equal(unname(g[i, w_top[i]]), 0)
            expect_gt(g[i, w_top[i] - 1], 0)
        }
    }
})

test_that("the upper boundary condition severs coupling above the growth cutoff", {
    # The transport coefficients sever the backward coupling from the
    # growth-chain top to the class above it (c = 0 there). Together with the
    # vanishing sub-diagonal (a = 0 at top+1, because growth out of the top
    # class is zero) this isolates the active spectrum from the inactive tail in
    # the tridiagonal solve, so the active solution never reads the tail. The
    # tail itself is then held at zero by zero_above_support(), which is what
    # stops a diffusive leak above the maximum size.
    p <- NS_params_small
    g <- getEGrowth(p)
    mu <- getMort(p)
    d <- getDiffusion(p)
    no_w <- length(p@w)
    w_top <- mizer:::support_top_idx(p)
    coefs <- mizer:::get_transport_coefs(p, p@initial_n, g, mu, dt = 1,
                                         recruitment_flux = getRDD(p), d = d,
                                         flux_limiter = "none")
    for (i in seq_len(nrow(p@species_params))) {
        # Backward sweep isolated: top class does not read the class above.
        expect_equal(unname(coefs$c[i, w_top[i]]), 0)
        # Forward sweep isolated: nothing grows out of the top class.
        if (w_top[i] < no_w) {
            expect_equal(unname(coefs$a[i, w_top[i] + 1]), 0)
        }
    }
})

test_that("the upper boundary condition stops diffusion leaking above w_max", {
    # On the full model (where the cutoff is a genuine maturity cutoff) a step
    # with diffusion leaves the density exactly zero above the growth-chain top.
    p <- setParams(NS_params, use_predation_diffusion = TRUE)
    no_w <- length(p@w)
    w_top <- mizer:::support_top_idx(p)
    sim <- project(p, t_max = 0.5, dt = 0.1, t_save = 0.5)
    nf <- finalN(sim)
    for (i in seq_len(nrow(p@species_params))) {
        if (w_top[i] < no_w) {
            expect_true(all(nf[i, (w_top[i] + 1):no_w] == 0))
        }
    }
})

test_that("steadyNewton handles initial guesses that are non-zero at large sizes where steady state is zero", {
    # Start with the stable model
    p <- p_steady
    # Fill the trailing tail (which is zero in p_steady) with positive numbers
    grid_top <- mizer:::support_top_idx(p)
    no_w <- length(p@w)

    # We will corrupt the tail of the first species with positive numbers
    # above its actual support.
    # In p_steady, species 1 (Sprat) only grows to its w_max (0.33g).
    # We set non-zero values up to the end of the grid.
    idx_zeros <- (grid_top[1] + 1):no_w
    if (length(idx_zeros) > 0) {
        p@initial_n[1, idx_zeros] <- 1e-3
    }

    # Run steadyNewton
    pn <- steadyNewton(p)
    expect_s4_class(pn, "MizerParams")

    # The tail should be correctly zeroed out in the result
    if (length(idx_zeros) > 0) {
        expect_true(all(pn@initial_n[1, idx_zeros] == 0))
    }
})
