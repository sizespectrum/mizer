# Tests for getLimitCycleSim()

# Use the same steadied model as test-steadyNewton.R
p_steady_lcs <- steady(NS_params_small, t_max = 50, progress_bar = FALSE)

test_that("getLimitCycleSim returns a MizerSim for a model with complex eigenvalues", {
    pn  <- steadyNewton(p_steady_lcs, stability = TRUE)
    stab <- attr(pn, "stability")

    # Only run test when the dominant eigenvalue is complex (Hopf mode dominant)
    skip_if(is.null(stab$hopf_period),
            "Model has no complex eigenvalues; skipping limit cycle test.")
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8,
            "Dominant eigenvalue is real; limit cycle test not applicable.")

    lcs <- getLimitCycleSim(pn)
    expect_s4_class(lcs, "MizerSim")
})

test_that("getLimitCycleSim time axis spans one period", {
    pn   <- steadyNewton(p_steady_lcs, stability = TRUE)
    stab <- attr(pn, "stability")

    skip_if(is.null(stab$hopf_period))
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8)

    lcs   <- getLimitCycleSim(pn)
    times <- getTimes(lcs)
    T_period <- stab$dominant_period

    expect_equal(times[1], 0)
    expect_lte(times[length(times)], ceiling(T_period) + 1e-8)
    expect_gt(length(times), 1L)
})

test_that("getLimitCycleSim abundances are non-negative", {
    pn   <- steadyNewton(p_steady_lcs, stability = TRUE)
    stab <- attr(pn, "stability")

    skip_if(is.null(stab$hopf_period))
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8)

    lcs <- getLimitCycleSim(pn, amplitude = 0.5)   # large amplitude stress test
    expect_true(all(lcs@n >= 0))
    expect_true(all(lcs@n_pp >= 0))
})

test_that("getLimitCycleSim t_save controls the time step spacing", {
    pn   <- steadyNewton(p_steady_lcs, stability = TRUE)
    stab <- attr(pn, "stability")

    skip_if(is.null(stab$hopf_period))
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8)

    T_period <- stab$dominant_period
    t_save   <- T_period / 50
    lcs      <- getLimitCycleSim(pn, t_save = t_save)
    times    <- getTimes(lcs)

    expect_equal(length(times), ceiling(T_period / t_save) + 1L)
    expect_equal(diff(times), rep(t_save, length(times) - 1L))
})

test_that("getLimitCycleSim n array has correct species and size dimnames", {
    pn   <- steadyNewton(p_steady_lcs, stability = TRUE)
    stab <- attr(pn, "stability")

    skip_if(is.null(stab$hopf_period))
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8)

    lcs <- getLimitCycleSim(pn)
    expect_equal(dimnames(lcs@n)$sp, dimnames(pn@initial_n)[[1]])
    expect_equal(dimnames(lcs@n)$w,  dimnames(pn@initial_n)[[2]])
})

test_that("getLimitCycleSim respects amplitude: max relative perturbation ~ amplitude", {
    pn   <- steadyNewton(p_steady_lcs, stability = TRUE)
    stab <- attr(pn, "stability")

    skip_if(is.null(stab$hopf_period))
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8)

    amp <- 0.1
    lcs <- getLimitCycleSim(pn, amplitude = amp, t_save = stab$dominant_period / 200)
    N_ss <- pn@initial_n
    active <- N_ss > 0
    max_rel <- max(abs(lcs@n - rep(N_ss, each = dim(lcs@n)[1])) /
                   rep(pmax(N_ss, .Machine$double.eps), each = dim(lcs@n)[1]),
                   na.rm = TRUE)
    # max relative perturbation should be close to amp (within floating-point)
    expect_equal(max_rel, amp, tolerance = 1e-4)
})
