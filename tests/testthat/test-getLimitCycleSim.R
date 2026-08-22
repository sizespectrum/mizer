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

    # A large amplitude drives cells negative, which is clipped and reported.
    expect_warning(lcs <- getLimitCycleSim(stab, amplitude = 0.5),
                   "clipped at zero")
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
    # n_pp(t) = n_pp* + A Re[e^{i omega t} v_resource], with A fixed by the
    # species biomass amplitude.
    bm <- function(n) as.numeric(sizeIntegral(pn, weighting = pn@w, n = n))
    B_ss <- bm(initialN(pn))
    c_sp <- complex(real = bm(Re(stab$hopf_eigenvector$fish)),
                    imaginary = bm(Im(stab$hopf_eigenvector$fish)))
    A <- lcs@sim_params$amplitude / max(Mod(c_sp[B_ss > 0]) / B_ss[B_ss > 0])
    omega <- Im(stab$hopf_eigenvalue)
    t2 <- getTimes(lcs)[2]
    expected <- initialNResource(pn) +
        A * Re(exp(1i * omega * t2) * stab$hopf_eigenvector$resource)
    expect_equal(as.numeric(npp[2, keep]), as.numeric(pmax(expected, 0)[keep]),
                 tolerance = 1e-10)
})

test_that("`amplitude` sets the largest relative swing in species biomass", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)
    skip_if(is.null(stab$hopf_eigenvalue))

    B_ss <- as.numeric(getBiomass(pn))
    for (amp in c(0.02, 0.05)) {
        lcs <- getLimitCycleSim(stab, amplitude = amp, t_save = 0.01)
        b   <- getBiomass(lcs)
        rel <- vapply(seq_along(B_ss),
                      function(i) max(abs(b[, i] - B_ss[i])) / B_ss[i],
                      numeric(1))
        # The hardest-swinging species reaches the cap; none exceeds it. The
        # time grid is discrete, so the peak is approached, not landed on.
        expect_lte(max(rel), amp * (1 + 1e-6))
        expect_gt(max(rel), amp * 0.99)
    }
})

test_that("`amplitude` scales the cycle linearly while nothing is clipped", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)
    skip_if(is.null(stab$hopf_eigenvalue))

    B_ss  <- as.numeric(getBiomass(pn))
    small <- getBiomass(getLimitCycleSim(stab, amplitude = 0.01))
    big   <- getBiomass(getLimitCycleSim(stab, amplitude = 0.02))
    dev_small <- sweep(small, 2, B_ss)
    dev_big   <- sweep(big, 2, B_ss)
    expect_equal(2 * as.numeric(dev_small), as.numeric(dev_big),
                 tolerance = 1e-8)
})

test_that("getLimitCycleSim reports clipping only when it matters", {
    skip_unless_experimental()
    pn   <- steadyNewton(p_steady_lcs)
    stab <- getStability(pn)
    skip_if(is.null(stab$hopf_eigenvalue))

    # A biomass amplitude does not bound the individual size classes, so
    # clipping sets in well below 1 and the report is the only thing that says
    # the picture has stopped being the linear mode.
    expect_silent(getLimitCycleSim(stab, amplitude = 0.001))
    expect_warning(getLimitCycleSim(stab, amplitude = 10),
                   "clipped at zero")
})
