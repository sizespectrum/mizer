# Tests for getOscillationModeSim()

# Every test here needs the same steady state and the same stability analysis,
# and that pair is by far the expensive part: an eigendecomposition of the whole
# coupled fish-resource Jacobian. They are built once, lazily, and shared. What
# is under test is getOscillationModeSim(), not the machinery that feeds it; the one
# test that deliberately starts from a raw MizerParams is marked as such.
delayedAssign("p_steady_lcs",
               suppressMessages(tuneSteadyState(NS_params_small, t_max = 50,
                                                progress_bar = FALSE)))
delayedAssign("lcs_params",
              suppressMessages(findSteadyState(p_steady_lcs,
                                               solver = "newton")))
delayedAssign("lcs_stab", suppressMessages(getStability(lcs_params)))

test_that("getOscillationModeSim returns a MizerSim for a model with complex eigenvalues", {
    skip_unless_experimental()
    pn  <- lcs_params
    stab <- lcs_stab

    # The function needs an oscillatory mode, which need not be the dominant
    # one: it draws `leading_oscillatory_eigenvector`, whose eigenvalue can sit well down the
    # spectrum below a real leading eigenvalue.
    skip_if(is.null(stab$leading_oscillatory_eigenvalue),
            "Model has no complex eigenvalues; skipping limit cycle test.")

    lcs <- getOscillationModeSim(stab)
    expect_s4_class(lcs, "MizerSim")
})

test_that("getOscillationModeSim time axis spans one period", {
    skip_unless_experimental()
    pn   <- lcs_params
    stab <- lcs_stab

    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    lcs   <- getOscillationModeSim(stab)
    times <- getTimes(lcs)
    T_period <- stab$oscillation_period

    expect_equal(times[1], 0)
    # Exactly at the period, not past it: that is what makes the last state the
    # first state again.
    expect_equal(times[length(times)], T_period)
    expect_gt(length(times), 1L)
})

test_that("getOscillationModeSim closes the cycle at a non-divisible period", {
    skip_unless_experimental()
    pn   <- lcs_params
    stab <- lcs_stab

    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    T_period <- stab$oscillation_period
    # Deliberately not a divisor of the period, which is the usual case: the
    # user picks a resolution, not a fraction of a period they do not yet know.
    t_save <- T_period / 7.3
    lcs    <- getOscillationModeSim(stab, t_save = t_save)
    times  <- getTimes(lcs)
    n_t    <- length(times)

    expect_equal(times[n_t], T_period)
    # Regular steps, with only the final interval shortened to land on T.
    expect_equal(diff(times)[-(n_t - 1)], rep(t_save, n_t - 2))
    expect_gt(diff(times)[n_t - 1], 0)
    expect_lte(diff(times)[n_t - 1], t_save * (1 + 1e-8))

    # The cycle closes: at t = T the phase factor is 1 again.
    expect_equal(lcs@n[n_t, , ], lcs@n[1, , ])
    expect_equal(lcs@n_pp[n_t, ], lcs@n_pp[1, ])
})

test_that("getOscillationModeSim validates its numerical controls", {
    skip_unless_experimental()
    stab <- lcs_stab
    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    expect_error(getOscillationModeSim(stab, amplitude = 0), "amplitude not greater")
    expect_error(getOscillationModeSim(stab, amplitude = -1), "amplitude not greater")
    expect_error(getOscillationModeSim(stab, amplitude = NA_real_), "amplitude")
    expect_error(getOscillationModeSim(stab, t_save = 0), "t_save not greater")
    expect_error(getOscillationModeSim(stab, t_save = Inf), "t_save")
})

test_that("getOscillationModeSim abundances are non-negative", {
    skip_unless_experimental()
    pn   <- lcs_params
    stab <- lcs_stab

    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    # A large amplitude drives cells negative, which is clipped and reported.
    expect_warning(lcs <- getOscillationModeSim(stab, amplitude = 0.5),
                   "clipped at zero")
    expect_true(all(lcs@n >= 0))
    expect_true(all(lcs@n_pp >= 0))
})

test_that("getOscillationModeSim t_save controls the time step spacing", {
    skip_unless_experimental()
    pn   <- lcs_params
    stab <- lcs_stab

    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    T_period <- stab$oscillation_period
    t_save   <- T_period / 50
    # The one test here that starts from a MizerParams rather than from the
    # cached stability list, so that path stays covered.
    lcs      <- getOscillationModeSim(pn, t_save = t_save)
    times    <- getTimes(lcs)

    expect_equal(times[length(times)], T_period)
    # Every interval but the last is exactly `t_save`; the last is whatever is
    # left over, which here is nothing to within rounding.
    expect_equal(diff(times)[-(length(times) - 1L)],
                 rep(t_save, length(times) - 2L))
    expect_lte(diff(times)[length(times) - 1L], t_save * (1 + 1e-8))
})

test_that("getOscillationModeSim n array has correct species and size dimnames", {
    skip_unless_experimental()
    pn   <- lcs_params
    stab <- lcs_stab

    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    lcs <- getOscillationModeSim(stab)
    expect_equal(dimnames(lcs@n)$sp, dimnames(pn@initial_n)[[1]])
    expect_equal(dimnames(lcs@n)$w,  dimnames(pn@initial_n)[[2]])
})

test_that("getOscillationModeSim oscillates the resource with the mode", {
    skip_unless_experimental()
    pn   <- lcs_params
    stab <- lcs_stab
    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    lcs <- getOscillationModeSim(stab, amplitude = 0.1)
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
    c_sp <- complex(real = bm(Re(stab$leading_oscillatory_eigenvector$fish)),
                    imaginary = bm(Im(stab$leading_oscillatory_eigenvector$fish)))
    A <- lcs@sim_params$amplitude / max(Mod(c_sp[B_ss > 0]) / B_ss[B_ss > 0])
    omega <- Im(stab$leading_oscillatory_eigenvalue)
    t2 <- getTimes(lcs)[2]
    expected <- initialNResource(pn) +
        A * Re(exp(1i * omega * t2) * stab$leading_oscillatory_eigenvector$resource)
    expect_equal(as.numeric(npp[2, keep]), as.numeric(pmax(expected, 0)[keep]),
                 tolerance = 1e-10)
})

test_that("`amplitude` sets the largest relative swing in species biomass", {
    skip_unless_experimental()
    pn   <- lcs_params
    stab <- lcs_stab
    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    B_ss <- as.numeric(getBiomass(pn))
    for (amp in c(0.02, 0.05)) {
        lcs <- getOscillationModeSim(stab, amplitude = amp, t_save = 0.01)
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
    pn   <- lcs_params
    stab <- lcs_stab
    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    B_ss  <- as.numeric(getBiomass(pn))
    small <- getBiomass(getOscillationModeSim(stab, amplitude = 0.01))
    big   <- getBiomass(getOscillationModeSim(stab, amplitude = 0.02))
    dev_small <- sweep(small, 2, B_ss)
    dev_big   <- sweep(big, 2, B_ss)
    expect_equal(2 * as.numeric(dev_small), as.numeric(dev_big),
                 tolerance = 1e-8)
})

test_that("getOscillationModeSim reports clipping only when it matters", {
    skip_unless_experimental()
    pn   <- lcs_params
    stab <- lcs_stab
    skip_if(is.null(stab$leading_oscillatory_eigenvalue))

    # A biomass amplitude does not bound the individual size classes, so
    # clipping sets in well below 1 and the report is the only thing that says
    # the picture has stopped being the linear mode.
    expect_silent(getOscillationModeSim(stab, amplitude = 0.001))
    expect_warning(getOscillationModeSim(stab, amplitude = 10),
                   "clipped at zero")
})
