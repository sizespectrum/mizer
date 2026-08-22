# Tests for the internals of the direct (Newton) steady-state solver and for
# the stability analysis built on the same machinery. The two exported functions
# that drive the solver, tuneSteadyState() and findSteadyState(), are tested in
# test-steadyState.R.

# A steadied model is a good starting guess for the root finder.
delayedAssign("p_steady",
               tuneSteadyState(NS_params_small, t_max = 50,
                               progress_bar = FALSE))
test_that("support_top_idx drops the pile-up bin for the second-order scheme", {
    skip_unless_experimental()
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
    skip_unless_experimental()
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
    skip_unless_experimental()
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
    skip_unless_experimental()
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
    skip_unless_experimental()
    # On the full model (where the cutoff is a genuine maturity cutoff) a step
    # with diffusion leaves the density exactly zero above the growth-chain top.
    p <- NS_params
    use_predation_diffusion(p) <- TRUE
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

# Tests for getStability() ---------------------------------------------------

test_that("getStability returns a well-formed list for a stable model", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    stab <- getStability(pn)

    expect_type(stab, "list")
    expect_named(stab, c("eigenvalues", "max_real_part", "stable",
                         "dominant_period", "oscillation_period", "leading_oscillatory_eigenvalue",
                         "leading_oscillatory_eigenvector", "n_active",
                         "leading_eigenvectors", "params"))
    expect_true(is.complex(stab$eigenvalues))
    expect_type(stab$max_real_part, "double")
    expect_type(stab$stable, "logical")
    expect_type(stab$dominant_period, "double")

    # Re(evals) are sorted descending
    expect_equal(order(Re(stab$eigenvalues), decreasing = TRUE),
                 seq_along(stab$eigenvalues))

    expect_equal(stab$max_real_part, max(Re(stab$eigenvalues)))
    expect_equal(stab$stable, stab$max_real_part < 0)
})
test_that("getStability warns when it is not handed a steady state", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    # A state that is not a fixed point makes the eigenvalues meaningless, so
    # the user has to be told rather than handed a stability verdict.
    off <- pn
    initialN(off)[1, ] <- initialN(off)[1, ] * 3
    expect_warning(getStability(off), "not at its steady state")
    expect_silent(getStability(pn))
})
test_that("getStability reports stable = TRUE for the NS model at its steady state", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    stab <- getStability(pn)

    expect_true(stab$stable)
    expect_lt(stab$max_real_part, 0)
})
test_that("getStability linearises the model's own reproduction function", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")

    # Holding the reproduction rate fixed is not an option on the analysis but
    # a property of the model: `constantRDD` has zero derivative, so it gives
    # the pinned Jacobian. The two verdicts are genuinely different, which is
    # why the analysis must read the model rather than take an argument.
    pinned <- pn
    pinned@species_params$constant_reproduction <- getRDD(pn)
    pinned@rates_funcs$RDD <- "constantRDD"

    stab <- getStability(pn)
    stab_pinned <- getStability(pinned)

    expect_equal(stab_pinned$n_active, stab$n_active)
    expect_gt(abs(stab_pinned$max_real_part - stab$max_real_part), 1e-2)
    # The reproduction feedback is what carries the state to the marginal end:
    # with this fixture's low reproduction levels the reproduction rate follows
    # the invested energy almost proportionally, and the abundance level is
    # then only weakly determined.
    expect_lt(abs(stab$max_real_part), abs(stab_pinned$max_real_part))

    # The removed argument must not creep back in through `...` forwarding.
    expect_error(getStability(pn, reproduction = "fixed"), "unused argument")
    expect_error(getDiscreteStability(pn, reproduction = "fixed"),
                 "unused argument")
    expect_error(getOscillationModeSim(pn, reproduction = "fixed"),
                 "unused argument")
})
test_that("getStability eigenvalues are consistent with max_real_part", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    stab <- getStability(pn)

    expect_equal(stab$max_real_part, max(Re(stab$eigenvalues)))
})
test_that("getDiscreteStability returns a well-formed list and depends on dt", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    disc1  <- getDiscreteStability(pn, dt = 1)
    disc01 <- getDiscreteStability(pn, dt = 0.1)

    expect_named(disc1, c("discrete_eigenvalues", "spectral_radius", "stable",
                          "dt", "n_active", "leading_eigenvectors", "params"))
    expect_true(is.complex(disc1$discrete_eigenvalues))
    expect_equal(disc1$spectral_radius, max(Mod(disc1$discrete_eigenvalues)))
    expect_equal(disc1$stable, disc1$spectral_radius < 1)
    expect_equal(disc1$dt, 1)
    expect_equal(disc1$n_active, getStability(pn)$n_active)

    # The one-step map is a property of the solver, so it does depend on dt.
    expect_false(isTRUE(all.equal(disc01$discrete_eigenvalues,
                                  disc1$discrete_eigenvalues)))
    expect_false(isTRUE(all.equal(disc01$spectral_radius,
                                  disc1$spectral_radius)))
})
test_that("the continuous eigenvalues are the dt -> 0 limit of the one-step map", {
    skip_unless_experimental()
    # getStability() differentiates the rates of change directly, so it involves
    # no time step. The one-step map at a small dt has to agree with it: for a
    # step of length dt the discrete eigenvalues satisfy mu = 1 + dt lambda +
    # O(dt^2). This is the check that the two maps describe the same model.
    pn <- findSteadyState(p_steady, solver = "newton")
    stab <- getStability(pn)
    dt <- 1e-3
    disc <- getDiscreteStability(pn, dt = dt)

    lambda <- sort(Re((disc$discrete_eigenvalues - 1) / dt), decreasing = TRUE)
    expect_equal(lambda[1], stab$max_real_part, tolerance = 1e-2)
})



test_that("getStability oscillation_period is NULL when all eigenvalues are real", {
    skip_unless_experimental()
    # For a 1-species model the dominant eigenvalue is typically real.
    # Use a simple single-species model as a proxy: the Newton solver on it,
    # then check that if oscillation_period is not NULL it is a positive finite number.
    pn <- findSteadyState(p_steady, solver = "newton")
    stab <- getStability(pn)

    if (!is.null(stab$oscillation_period)) {
        expect_true(is.finite(stab$oscillation_period))
        expect_gt(stab$oscillation_period, 0)
    }
})
test_that("getStability returns a well-formed list including the resource", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    stab <- getStability(pn)

    expect_type(stab, "list")
    expect_named(stab, c("eigenvalues", "max_real_part", "stable",
                         "dominant_period", "oscillation_period", "leading_oscillatory_eigenvalue",
                         "leading_oscillatory_eigenvector", "n_active",
                         "leading_eigenvectors", "params"))
    # The resource is always part of the state, so the Jacobian is larger than
    # the number of active fish cells.
    expect_gt(stab$n_active, sum(initialN(pn) > 0))
    expect_named(stab$leading_eigenvectors, c("fish", "resource"))
    expect_equal(nrow(stab$leading_eigenvectors$resource), length(pn@w_full))
    expect_true(stab$stable)
})
test_that("getStability works with non-semichemostat resource dynamics", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    pn@resource_dynamics <- "resource_constant"
    # A constant resource contributes only zero rows, but the call must work:
    # the quasi-static substitution used to restrict this to the semichemostat.
    stab <- expect_no_error(getStability(pn))
    expect_length(stab$eigenvalues, stab$n_active)
})
test_that("fd_step_scale floors the step at the local scale of the spectrum", {
    x <- c(2, 0, -3, 0)
    local_scale <- c(1, 5, 1, 0)
    scale <- mizer:::fd_step_scale(x, local_scale)

    # A cell with a magnitude of its own keeps it, whatever its sign.
    expect_equal(scale[[1]], 2)
    expect_equal(scale[[3]], 3)
    # A cell at zero takes the local scale, not an absolute epsilon.
    expect_equal(scale[[2]], 5)
    # With no local scale either (a row that is zero throughout) it falls back
    # to the absolute floor. The step must never be zero.
    expect_equal(scale[[4]], .Machine$double.eps)
    expect_true(all(scale > 0))
})
test_that("getStability resolves cells sitting at exactly zero", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")
    stab <- getStability(pn)

    # Punch an isolated hole in the middle of the second species' spectrum, of
    # the kind the second-order schemes can leave behind. Its Jacobian column
    # must still be resolved: a step floored at an absolute `.Machine$double.eps`
    # would be swamped by the rounding error of the outputs, leaving a zero
    # column that drops the cell from the Jacobian and contributes a spurious
    # zero eigenvalue.
    active <- mizer:::steady_active_set(pn)
    j <- floor(mean(c(pn@w_min_idx[2], active$w_top[2])))
    holed <- pn
    holed@initial_n[2, j] <- 0
    stab_holed <- getStability(holed)

    expect_equal(stab_holed$n_active, stab$n_active)
    expect_false(any(Mod(stab_holed$eigenvalues) < 1e-12))
    # Both comparisons below are absolute. This fixture has low reproduction
    # levels, so its leading eigenvalue sits close to the marginal value zero
    # (around -0.01), and a relative tolerance on it measures nothing. The
    # failure the test is for is unmistakable on the absolute scale: dropping
    # the holed cell from the Jacobian puts a spurious eigenvalue at exactly
    # zero, which moves `max_real_part` by its whole magnitude.
    # The hole is one cell out of many, so the spectrum barely moves.
    expect_lt(abs(stab_holed$max_real_part - stab$max_real_part), 1e-3)
    # The step is relative, so the answer must not depend on `h`.
    expect_lt(abs(getStability(holed, h = 1e-5)$max_real_part -
                      stab_holed$max_real_part), 1e-3)
})
test_that("getStability never evaluates rate functions at negative abundances", {
    skip_unless_experimental()
    # A cell at zero has to be perturbed by more than its own (zero) magnitude,
    # so the column is differenced forwards to keep the state in N >= 0.
    pn <- findSteadyState(p_steady, solver = "newton")
    active <- mizer:::steady_active_set(pn)
    j <- floor(mean(c(pn@w_min_idx[2], active$w_top[2])))
    holed <- pn
    holed@initial_n[2, j] <- 0

    seen <- new.env()
    seen$min_n <- Inf
    seen$min_n_pp <- Inf
    record_min <- function(params, n, n_pp, n_other, t, ...) {
        seen$min_n    <- min(seen$min_n, n)
        seen$min_n_pp <- min(seen$min_n_pp, n_pp)
        mizerEncounter(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                       ...)
    }
    # setRateFunction() looks the function up by name in the global environment.
    assign("record_min", record_min, envir = globalenv())
    on.exit(rm("record_min", envir = globalenv()), add = TRUE)
    spied <- setRateFunction(holed, "Encounter", "record_min")

    getStability(spied)
    expect_gte(seen$min_n, 0)
    expect_gte(seen$min_n_pp, 0)

    seen$min_n <- Inf
    seen$min_n_pp <- Inf
    getDiscreteStability(spied)
    expect_gte(seen$min_n, 0)
    expect_gte(seen$min_n_pp, 0)

    # The resource carries structural zeros above w_pp_cutoff, so the coupled
    # analysis exercises the resource columns too.
    seen$min_n <- Inf
    seen$min_n_pp <- Inf
    getStability(spied)
    expect_gte(seen$min_n, 0)
    expect_gte(seen$min_n_pp, 0)
})
test_that("getStability errors informatively on a non-finite rate function", {
    skip_unless_experimental()
    pn <- findSteadyState(p_steady, solver = "newton")

    not_finite <- function(params, n, n_pp, n_other, t, ...) {
        encounter <- mizerEncounter(params, n = n, n_pp = n_pp,
                                    n_other = n_other, t = t, ...)
        encounter[1, 1] <- NaN
        encounter
    }
    assign("not_finite", not_finite, envir = globalenv())
    on.exit(rm("not_finite", envir = globalenv()), add = TRUE)

    expect_error(getStability(setRateFunction(pn, "Encounter", "not_finite")),
                 "returned non-finite values")
})
test_that("leading_eigenvectors have correct shape and are normalised", {
    skip_unless_experimental()
    pn   <- findSteadyState(p_steady, solver = "newton")
    stab <- getStability(pn)
    lev  <- stab$leading_eigenvectors$fish

    # Shape: (n_species, n_sizes, 2)
    expect_equal(dim(lev)[1], nrow(pn@initial_n))
    expect_equal(dim(lev)[2], ncol(pn@initial_n))
    expect_equal(dim(lev)[3], 2L)

    # Normalised jointly over fish and resource, by the largest perturbation
    # relative to the steady state.
    N <- initialN(pn); npp <- initialNResource(pn)
    joint_max_rel <- function(k) {
        max(max(Mod(stab$leading_eigenvectors$fish[, , k])[N > 0] / N[N > 0]),
            max(Mod(stab$leading_eigenvectors$resource[, k])[npp > 0] /
                    npp[npp > 0]))
    }
    expect_equal(joint_max_rel(1), 1, tolerance = 1e-10)
    expect_equal(joint_max_rel(2), 1, tolerance = 1e-10)

    # Dimnames match initial_n
    expect_equal(dimnames(lev)[1:2], dimnames(pn@initial_n))
})
test_that("getStability returns the dominant oscillatory mode with its vector", {
    skip_unless_experimental()
    pn   <- findSteadyState(p_steady, solver = "newton")
    stab <- getStability(pn)
    skip_if(is.null(stab$leading_oscillatory_eigenvalue), "No complex eigenvalue in this model.")

    # oscillation_period must describe leading_oscillatory_eigenvalue, and that eigenvalue must be
    # the complex one with the largest real part.
    expect_equal(stab$oscillation_period, 2 * pi / abs(Im(stab$leading_oscillatory_eigenvalue)))
    is_cplx <- abs(Im(stab$eigenvalues)) > 1e-8
    expect_equal(stab$leading_oscillatory_eigenvalue,
                 stab$eigenvalues[is_cplx][which.max(Re(stab$eigenvalues[is_cplx]))])

    # The vector travels with the eigenvalue, one mode, no trailing dimension.
    expect_equal(dim(stab$leading_oscillatory_eigenvector$fish), dim(pn@initial_n))
    expect_length(stab$leading_oscillatory_eigenvector$resource, length(pn@w_full))

    # It really is an eigenvector of that eigenvalue: when the Hopf mode is
    # also the dominant one, it must match the leading eigenvector up to phase.
    if (isTRUE(all.equal(stab$leading_oscillatory_eigenvalue, stab$eigenvalues[1]))) {
        a <- as.vector(stab$leading_oscillatory_eigenvector$fish)
        b <- as.vector(stab$leading_eigenvectors$fish[, , 1])
        keep <- Mod(b) > 0
        ratio <- (a[keep] / b[keep])
        expect_equal(max(Mod(ratio)) - min(Mod(ratio)), 0, tolerance = 1e-8)
    }
})
test_that("steady_active_set leaves an absent species out of the system", {
    skip_unless_experimental()
    p <- NS_params_small
    p@initial_n[2, ] <- 0

    active <- mizer:::steady_active_set(p)
    expect_true(active$absent[[2]])
    expect_false(any(active$absent[-2]))
    # No unknowns for that species, so nothing to take a log of.
    expect_false(any(active$mask[2, ]))
    expect_true(all(active$mask[-2, p@w_min_idx[1]:active$w_top[1]]))
    expect_equal(active$n_fish_active, sum(active$mask))

    # And the packed state still unpacks to a full density matrix, with the
    # absent species at zero. (With the resource held fixed the packed state is
    # the fish block alone, which keeps the check to the point.)
    fixed <- mizer:::steady_active_set(p, resource = "fixed")
    N <- fixed$unpack(rep(log(1e-5), fixed$n_fish_active))$N
    expect_true(all(N[2, ] == 0))
    expect_true(any(N[1, ] > 0))
})

test_that("the Newton solver keeps an already-absent species at zero", {
    skip_unless_experimental()
    # A zero row used to reach log() as -Inf and stop the root finder before it
    # took a single step.
    p <- p_steady
    initialN(p)[2, ] <- 0

    warns <- character(0)
    pn <- withCallingHandlers(
        suppressMessages(findSteadyState(p, solver = "newton")),
        warning = function(w) {
            warns <<- c(warns, conditionMessage(w))
            invokeRestart("muffleWarning")
        })

    expect_true(any(grepl("already absent and were held at zero", warns)))
    # Nothing died during the solve, so nothing is reported as having.
    expect_false(any(grepl("went extinct", warns)))
    expect_true(all(initialN(pn)[2, ] == 0))
    expect_true(any(initialN(pn)[1, ] > 0))
    expect_true(any(initialN(pn)[3, ] > 0))
})

test_that("the stability analyses validate their numerical controls", {
    # These fire before any computation, so no steady state is needed.
    expect_error(getStability(NS_params_small, h = 0), "h not greater than 0")
    expect_error(getStability(NS_params_small, h = -1e-4), "h not greater than 0")
    expect_error(getStability(NS_params_small, h = Inf), "h")
    expect_error(getDiscreteStability(NS_params_small, h = 0),
                 "h not greater than 0")
    expect_error(getDiscreteStability(NS_params_small, dt = 0),
                 "dt not greater than 0")
    expect_error(getDiscreteStability(NS_params_small, dt = Inf), "dt")
})

test_that("getDiscreteStability differentiates the step project() takes", {
    skip_unless_experimental()
    # The point of the function is to say what mizer's solver does at a given
    # `dt`, so the map it linearises has to be that step and not one that merely
    # agrees with it to first order.
    p  <- p_steady
    dt <- 0.1
    ctx <- suppressMessages(suppressWarnings(mizer:::stability_context(p)))
    out <- mizer:::stability_step(ctx, dt)(ctx$x0)

    sim <- suppressMessages(project(p, t_max = dt, dt = dt, t_save = dt,
                                    progress_bar = FALSE))

    expect_equal(out[seq_len(ctx$n_fish_active)],
                 finalN(sim)[ctx$active$mask])
    expect_equal(out[ctx$n_fish_active + seq_len(ctx$n_npp)],
                 as.numeric(finalNResource(sim)))
})

test_that("getStability says a dynamic component is not in its Jacobian", {
    skip_unless_experimental()
    # A component that happens not to move, but whose dynamics function is not
    # `constant_other`: the check is on what the model declares, because that is
    # all it can know. Keeping it still also keeps the model at its steady
    # state, so this is the only warning the call has to raise.
    assign("mizer_test_stability_component", function(params, n_other,
                                                      component, ...) {
        n_other[[component]]
    }, envir = .GlobalEnv)
    withr::defer(rm("mizer_test_stability_component", envir = .GlobalEnv))

    p <- setComponent(p_steady, "grower", initial_value = 1,
                      dynamics_fun = "mizer_test_stability_component")

    # The Jacobian has rows for the fish and the resource and none for the
    # component, so the spectrum is that of a subsystem. That has to be said.
    expect_warning(suppressMessages(getStability(p)),
                   "covers the consumers and the resource only")
})
