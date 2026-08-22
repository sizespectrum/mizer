# Tests for getLimitCycleSim()

# Use the same steadied model as test-steadyNewton.R
delayedAssign("p_steady_lcs",
               steady(NS_params_small, t_max = 50, progress_bar = FALSE))

test_that("getLimitCycleSim returns a MizerSim for a model with complex eigenvalues", {
    skip_unless_experimental()
    pn  <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)

    # Only run test when the dominant eigenvalue is complex (Hopf mode dominant)
    skip_if(is.null(stab$hopf_period),
            "Model has no complex eigenvalues; skipping limit cycle test.")
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8,
            "Dominant eigenvalue is real; limit cycle test not applicable.")

    lcs <- getLimitCycleSim(stab)
    expect_s4_class(lcs, "MizerSim")
})

test_that("getLimitCycleSim time axis spans one period", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)

    skip_if(is.null(stab$hopf_period))
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8)

    lcs   <- getLimitCycleSim(stab)
    times <- getTimes(lcs)
    T_period <- stab$dominant_period

    expect_equal(times[1], 0)
    expect_lte(times[length(times)], ceiling(T_period) + 1e-8)
    expect_gt(length(times), 1L)
})

test_that("getLimitCycleSim abundances are non-negative", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)

    skip_if(is.null(stab$hopf_period))
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8)

    lcs <- getLimitCycleSim(stab, amplitude = 0.5)   # large amplitude stress test
    expect_true(all(lcs@n >= 0))
    expect_true(all(lcs@n_pp >= 0))
})

test_that("getLimitCycleSim t_save controls the time step spacing", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)

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
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)

    skip_if(is.null(stab$hopf_period))
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8)

    lcs <- getLimitCycleSim(stab)
    expect_equal(dimnames(lcs@n)$sp, dimnames(pn@initial_n)[[1]])
    expect_equal(dimnames(lcs@n)$w,  dimnames(pn@initial_n)[[2]])
})

test_that("getLimitCycleSim respects amplitude: max relative perturbation ~ amplitude", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)

    skip_if(is.null(stab$hopf_period))
    skip_if(abs(Im(stab$eigenvalues[1])) <= 1e-8)

    amp <- 0.1
    lcs <- getLimitCycleSim(stab, amplitude = amp, t_save = stab$dominant_period / 200)
    N_ss <- pn@initial_n
    active <- N_ss > 0
    max_rel <- max(abs(lcs@n - rep(N_ss, each = dim(lcs@n)[1])) /
                   rep(pmax(N_ss, .Machine$double.eps), each = dim(lcs@n)[1]),
                   na.rm = TRUE)
    # max relative perturbation should be close to amp (within floating-point)
    expect_equal(max_rel, amp, tolerance = 1e-4)
})

test_that("getLimitCycleSim oscillates the resource with the mode", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)
    skip_if(is.null(stab$hopf_eigenvalue))

    lcs <- getLimitCycleSim(stab, amplitude = 0.1)
    npp <- NResource(lcs)
    keep <- npp[1, ] > 0

    # The resource must actually move: it used to be slaved to the fish at
    # their quasi-static equilibrium, which threw away its phase.
    swing <- apply(npp[, keep, drop = FALSE], 2,
                   function(x) (max(x) - min(x)) / mean(x))
    expect_gt(max(swing), 0)

    # And it must move as the eigenvector says, not as a function of the fish:
    # n_pp(t) = n_pp* + A Re[e^{i omega t} v_resource].
    A <- lcs@sim_params$amplitude / max(
        max(Mod(stab$hopf_eigenvector$fish)[initialN(pn) > 0] /
                initialN(pn)[initialN(pn) > 0]),
        max(Mod(stab$hopf_eigenvector$resource)[keep] /
                initialNResource(pn)[keep]))
    omega <- Im(stab$hopf_eigenvalue)
    t2 <- getTimes(lcs)[2]
    expected <- initialNResource(pn) +
        A * Re(exp(1i * omega * t2) * stab$hopf_eigenvector$resource)
    expect_equal(as.numeric(npp[2, keep]), as.numeric(pmax(expected, 0)[keep]),
                 tolerance = 1e-10)
})

test_that("getLimitCycleSim caps the relative perturbation at `amplitude`", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)
    skip_if(is.null(stab$hopf_eigenvalue))

    for (amp in c(0.05, 0.2)) {
        lcs <- getLimitCycleSim(stab, amplitude = amp)
        N   <- N(lcs)
        N_ss <- initialN(pn)
        rel <- apply(N, 1, function(m) {
            max(abs(m[N_ss > 0] - N_ss[N_ss > 0]) / N_ss[N_ss > 0])
        })
        # The cap is over the whole state, so the fish reach it only when the
        # mode puts its largest relative swing on them; never exceed it.
        expect_lte(max(rel), amp + 1e-8)
    }
})

test_that("getLimitCycleSim reports clipping only when it matters", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)
    skip_if(is.null(stab$hopf_eigenvalue))

    # Below 1 the perturbation is smaller than the state it perturbs, so no
    # positive cell can go negative and there is nothing to report.
    expect_silent(getLimitCycleSim(stab, amplitude = 0.1))
    # Above 1 it can, and must be.
    expect_warning(getLimitCycleSim(stab, amplitude = 3),
                   "clipped at zero")
})
