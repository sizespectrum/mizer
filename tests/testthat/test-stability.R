delayedAssign("p_stability",
              findSteadyState(tuneSteadyState(
                  NS_params_small, t_max = 50, progress_bar = FALSE
              ), solver = "newton"))

test_that("getStability returns a well-formed list for a stable model", {
    skip_unless_experimental()
    pn <- p_stability
    stab <- getStability(pn)

    expect_type(stab, "list")
    expect_named(stab, c(
        "eigenvalues", "max_real_part", "dominant_period", "oscillation_period",
        "leading_oscillatory_eigenvalue", "leading_oscillatory_eigenvector",
        "n_active", "leading_eigenvectors", "params", "status", "stable",
        "scope", "fixed", "prescribed", "prescribed_masks", "coverage",
        "autonomous", "time", "tolerance", "assessment"
    ))
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
test_that("non-equilibrium states have diagnostics without a verdict", {
    skip_unless_experimental()
    pn <- p_stability
    # A state that is not a fixed point makes the eigenvalues meaningless, so
    # the user has to be told rather than handed a stability verdict.
    off <- pn
    initialN(off)[1, ] <- initialN(off)[1, ] * 3
    result <- getStability(off)
    expect_false(result$assessment$verdict)
    expect_identical(result$status, "inconclusive")
    expect_true(is.na(result$stable))
    expect_true(is.complex(result$eigenvalues))
    expect_silent(getStability(pn))
})
test_that("the NS equilibrium is strictly stable", {
    skip_unless_experimental()
    pn <- p_stability
    stab <- getStability(pn)

    expect_true(stab$stable)
    expect_lt(stab$max_real_part, 0)
})
test_that("getStability linearises the model's own reproduction function", {
    skip_unless_experimental()
    pn <- p_stability

    # Holding the reproduction rate fixed is not an option on the analysis but
    # a property of the model: `constantRDD` has zero derivative, so it gives
    # the pinned Jacobian. The two verdicts are genuinely different, which is
    # why the analysis must read the model rather than take an argument.
    pinned <- pn
    pinned@species_params$constant_reproduction <- getRDD(pn)
    pinned <- setRateFunction(pinned, "RDD", "constantRDD")

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
    pn <- p_stability
    stab <- getStability(pn)

    expect_equal(stab$max_real_part, max(Re(stab$eigenvalues)))
})
test_that("getDiscreteStability returns a well-formed list and depends on dt", {
    skip_unless_experimental()
    pn <- p_stability
    disc1  <- getDiscreteStability(pn, dt = 1)
    disc01 <- getDiscreteStability(pn, dt = 0.1)

    expect_named(disc1, c(
        "discrete_eigenvalues", "spectral_radius", "dt", "n_active",
        "leading_eigenvectors", "params", "status", "stable", "scope", "fixed",
        "prescribed", "prescribed_masks", "coverage", "autonomous", "time",
        "tolerance", "assessment"
    ))
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
test_that("discrete eigenvalues converge to the continuous limit", {
    skip_unless_experimental()
    # getStability() differentiates the rates of change directly, so it involves
    # no time step. The one-step map at a small dt has to agree with it: for a
    # step of length dt the discrete eigenvalues satisfy mu = 1 + dt lambda +
    # O(dt^2). This is the check that the two maps describe the same model.
    pn <- p_stability
    stab <- getStability(pn)
    dt <- 1e-3
    disc <- getDiscreteStability(pn, dt = dt)

    lambda <- sort(Re((disc$discrete_eigenvalues - 1) / dt), decreasing = TRUE)
    expect_equal(lambda[1], stab$max_real_part, tolerance = 1e-2)
})



test_that("reported oscillation periods are positive and finite", {
    skip_unless_experimental()
    # For a 1-species model the dominant eigenvalue is typically real.
    # Use a simple single-species model as a proxy: the Newton solver on it,
    # Any reported oscillation_period must be positive and finite.
    pn <- p_stability
    stab <- getStability(pn)

    if (!is.null(stab$oscillation_period)) {
        expect_true(is.finite(stab$oscillation_period))
        expect_gt(stab$oscillation_period, 0)
    }
})
test_that("getStability returns a well-formed list including the resource", {
    skip_unless_experimental()
    pn <- p_stability
    stab <- getStability(pn)

    expect_type(stab, "list")
    expect_named(stab, c(
        "eigenvalues", "max_real_part", "dominant_period", "oscillation_period",
        "leading_oscillatory_eigenvalue", "leading_oscillatory_eigenvector",
        "n_active", "leading_eigenvectors", "params", "status", "stable",
        "scope", "fixed", "prescribed", "prescribed_masks", "coverage",
        "autonomous", "time", "tolerance", "assessment"
    ))
    # The resource is always part of the state, so the Jacobian is larger than
    # the number of active fish cells.
    expect_gt(stab$n_active, sum(initialN(pn) > 0))
    expect_named(stab$leading_eigenvectors, c("fish", "resource"))
    expect_equal(nrow(stab$leading_eigenvectors$resource), length(pn@w_full))
    expect_true(stab$stable)
})
test_that("getStability works with non-semichemostat resource dynamics", {
    skip_unless_experimental()
    pn <- p_stability
    pn <- setResource(pn, resource_dynamics = "resource_constant",
                      balance = FALSE)
    # A known constant resource is an explicitly prescribed input, so its
    # coordinates are excluded and its reconstructed eigenvector entries zero.
    stab <- expect_no_error(getStability(pn))
    expect_length(stab$eigenvalues, stab$n_active)
    expect_identical(stab$prescribed, "resource")
    expect_true(all(stab$leading_eigenvectors$resource == 0))
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
    pn <- p_stability
    stab <- getStability(pn)

    # Punch an isolated hole in the middle of the second species' spectrum, of
    # the kind the second-order schemes can leave behind. Its Jacobian column
    # must still be resolved: an absolute `.Machine$double.eps` step
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
test_that("stability perturbations have nonnegative abundances", {
    skip_unless_experimental()
    # A cell at zero has to be perturbed by more than its own (zero) magnitude,
    # so the column is differenced forwards to keep the state in N >= 0.
    pn <- p_stability
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
    spied <- setProcess(spied, "consumers", autonomous = TRUE)
    spied <- setProcess(spied, "resource", autonomous = TRUE)

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
    pn <- p_stability

    not_finite <- function(params, n, n_pp, n_other, t, ...) {
        encounter <- mizerEncounter(params, n = n, n_pp = n_pp,
                                    n_other = n_other, t = t, ...)
        encounter[1, 1] <- NaN
        encounter
    }
    assign("not_finite", not_finite, envir = globalenv())
    on.exit(rm("not_finite", envir = globalenv()), add = TRUE)

    pn@rates_funcs$Encounter <- "not_finite"
    pn <- setProcess(pn, "consumers", autonomous = TRUE)
    pn <- setProcess(pn, "resource", autonomous = TRUE)
    expect_error(getStability(pn), "finite")
})
test_that("leading_eigenvectors have correct shape and are normalised", {
    skip_unless_experimental()
    pn   <- p_stability
    stab <- getStability(pn)
    lev  <- stab$leading_eigenvectors$fish

    # Shape: (n_species, n_sizes, 2)
    expect_equal(dim(lev)[1], nrow(pn@initial_n))
    expect_equal(dim(lev)[2], ncol(pn@initial_n))
    expect_equal(dim(lev)[3], 2L)

    # Normalised jointly over fish and resource, by the largest perturbation
    # relative to the steady state.
    n <- initialN(pn)
    npp <- initialNResource(pn)
    joint_max_rel <- function(k) {
        max(max(Mod(stab$leading_eigenvectors$fish[, , k])[n > 0] / n[n > 0]),
            max(Mod(stab$leading_eigenvectors$resource[, k])[npp > 0] /
                    npp[npp > 0]))
    }
    expect_equal(joint_max_rel(1), 1, tolerance = 1e-10)
    expect_equal(joint_max_rel(2), 1, tolerance = 1e-10)

    # Dimnames match initial_n
    expect_equal(dimnames(lev)[1:2], dimnames(pn@initial_n))
})
test_that("the dominant oscillatory mode carries its vector", {
    skip_unless_experimental()
    pn   <- p_stability
    stab <- getStability(pn)
    skip_if(is.null(stab$leading_oscillatory_eigenvalue),
            "No complex eigenvalue in this model.")

    # The period belongs to the oscillatory eigenvalue with largest real part.
    value <- stab$leading_oscillatory_eigenvalue
    vector <- stab$leading_oscillatory_eigenvector
    expect_equal(stab$oscillation_period, 2 * pi / abs(Im(value)))
    is_cplx <- abs(Im(stab$eigenvalues)) > 1e-8
    complex_values <- stab$eigenvalues[is_cplx]
    expect_equal(value, complex_values[which.max(Re(complex_values))])

    # The vector travels with the eigenvalue, one mode, no trailing dimension.
    expect_equal(dim(vector$fish), dim(pn@initial_n))
    expect_length(vector$resource, length(pn@w_full))

    # It really is an eigenvector of that eigenvalue: when the Hopf mode is
    # also the dominant one, it must match the leading eigenvector up to phase.
    if (isTRUE(all.equal(value, stab$eigenvalues[1]))) {
        a <- as.vector(stab$leading_oscillatory_eigenvector$fish)
        b <- as.vector(stab$leading_eigenvectors$fish[, , 1])
        keep <- Mod(b) > 0
        ratio <- (a[keep] / b[keep])
        expect_equal(max(Mod(ratio)) - min(Mod(ratio)), 0, tolerance = 1e-8)
    }
})
test_that("the stability analyses validate their numerical controls", {
    # These fire before any computation, so no steady state is needed.
    expect_error(getStability(NS_params_small, h = 0), "h not greater than 0")
    expect_error(getStability(NS_params_small, h = -1e-4),
                 "h not greater than 0")
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
    p  <- p_stability
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
    # all it can know. Without a convergence declaration its steady-state
    # assessment is also incomplete, even though the callback returns its input.
    assign("mizer_test_stability_component", function(params, n_other,
                                                      component, ...) {
        n_other[[component]]
    }, envir = .GlobalEnv)
    withr::defer(rm("mizer_test_stability_component", envir = .GlobalEnv))

    p <- setComponent(p_stability, "grower", initial_value = 1,
                      dynamics_fun = "mizer_test_stability_component")

    whole <- getStability(p)
    expect_identical(whole$status, "inconclusive")
    expect_true(is.na(whole$stable))
    expect_null(whole$eigenvalues)
    expect_true(is.na(whole$n_active))
    expect_contains(whole$coverage$missing, "component:grower")
    conditional <- getStability(p, scope = "conditional",
                                fixed = "component:grower")
    expect_true(conditional$coverage$complete)
    expect_true(conditional$stable)
    expect_identical(conditional$fixed, "component:grower")
    expect_identical(conditional$scope, "conditional")
})

test_that("stability maps give custom resources the complete effective rates", {
    calls <- 0L
    rlang::local_bindings(
        stability_complete_rates = function(params, ...) {
            calls <<- calls + 1L
            rates <- mizerRates(params, ...)
            rates$rdd <- rates$rdd / 2
            rates
        },
        stability_fishing_equation = function(params, n_pp, rates, t) {
            n_pp * sum(rates$f_mort)
        },
        stability_fishing_resource = function(n_pp, rates, dt, predicted_state,
                                              endpoint_rates, ...) {
            if (is.null(rates$f_mort)) stop("missing fishing mortality")
            stopifnot(is.null(predicted_state), is.null(endpoint_rates))
            n_pp * (1 + dt * sum(rates$f_mort))
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    p@resource_dynamics <- "stability_fishing_resource"
    p@rates_funcs$Rates <- "stability_complete_rates"
    p <- setProcess(p, "consumers", autonomous = TRUE)
    p <- setProcess(p, "resource", role = "dynamic", autonomous = TRUE,
                    rate_fun = "stability_fishing_equation")
    ctx <- suppressMessages(suppressWarnings(stability_context(p)))
    rates <- mizerRates(p, n = p@initial_n, n_pp = p@initial_n_pp,
                        n_other = p@initial_n_other, effort = p@initial_effort,
                        rates_fns = projectRateFunctions(p))
    fishing <- sum(rates$f_mort)
    expect_gt(fishing, 0)

    calls <- 0L
    rhs <- stability_rhs(ctx)(ctx$x0)
    expect_identical(calls, 1L)
    expected <- -consumer_residual(p, n = p@initial_n,
                                   n_pp = p@initial_n_pp,
                                   n_other = p@initial_n_other,
                                   effort = p@initial_effort,
                                   rates = rates, rdd = rates$rdd / 2)
    expect_equal(rhs[seq_len(ctx$n_fish_active)], expected[ctx$active$mask])
    resource_idx <- ctx$n_fish_active + seq_len(ctx$n_npp)
    present <- p@initial_n_pp > 0
    expect_equal(rhs[resource_idx[present]] /
                     as.numeric(p@initial_n_pp[present]),
                 rep(fishing, sum(present)), tolerance = 1e-10)
    expect_true(all(rhs[resource_idx[!present]] == 0))

    calls <- 0L
    dt <- 0.001
    stepped <- stability_step(ctx, dt)(ctx$x0)
    expect_identical(calls, 1L)
    # Compare the map with the public Euler projection, including the custom
    # provider's assembled recruitment and the same resource update.
    sim <- project(p, t_max = dt, dt = dt, t_save = dt, method = "euler",
                   progress_bar = FALSE, info_level = 0)
    expect_equal(stepped[seq_len(ctx$n_fish_active)],
                 finalN(sim)[ctx$active$mask])
    expect_equal(stepped[resource_idx], as.numeric(finalNResource(sim)))
})


test_that("absent species retain invasion coordinates in both Jacobians", {
    rlang::local_bindings(
        stability_invasion_rates = function(params, n, ...) {
            rates <- mizerRates(params, n = n, ...)
            rates$e_growth[] <- 0
            rates$diffusion[] <- 0
            rates$mort[] <- 0.5
            rates$rdd <- vapply(seq_len(nrow(n)), function(i) {
                j <- params$w_min_idx[i]
                2 * n[i, j] * params$dw[j]
            }, numeric(1))
            rates
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    p$initial_n[] <- 0
    p <- setResource(p, resource_dynamics = "resource_constant",
                     balance = FALSE)
    p <- setRateFunction(p, "Rates", "stability_invasion_rates")
    p <- setProcess(p, "consumers", autonomous = TRUE)
    p <- setProcess(p, "resource", autonomous = TRUE)
    expect_false(any(steady_active_set(p)$mask))
    continuous <- getStability(p)
    discrete <- getDiscreteStability(p, dt = 0.1)
    expected <- sum(support_top_idx(p) - p$w_min_idx + 1L)
    expect_equal(continuous$n_active, expected)
    expect_equal(discrete$n_active, expected)
    expect_equal(continuous$max_real_part, 1.5, tolerance = 1e-8)
    expect_equal(discrete$spectral_radius, 1.2 / 1.05, tolerance = 1e-8)
    expect_identical(continuous$status, "unstable")
    expect_identical(discrete$status, "unstable")
    expect_false(continuous$stable)
})

test_that("zero turnover and capacity retain neutral resource modes", {
    p <- NS_params_small
    p$initial_n[] <- 0
    p$initial_n_pp[] <- 0
    p$rr_pp[] <- 0
    p$cc_pp[] <- 0
    continuous <- getStability(p, scope = "conditional", fixed = "consumers")
    discrete <- getDiscreteStability(p, scope = "conditional",
                                     fixed = "consumers")
    expect_true(continuous$assessment$verdict)
    expect_true(discrete$assessment$verdict)
    expect_equal(continuous$n_active, length(p$w_full))
    expect_equal(Re(continuous$eigenvalues), rep(0, length(p$w_full)))
    expect_equal(Re(discrete$discrete_eigenvalues), rep(1, length(p$w_full)))
    expect_identical(continuous$status, "inconclusive")
    expect_identical(discrete$status, "inconclusive")
    expect_true(is.na(continuous$stable))
    expect_true(is.na(discrete$stable))
})

test_that("partial prescribed masks precede linearisation and preserve zeros", {
    rlang::local_bindings(
        stability_resource_input = function(params) {
            mask <- params$initial_n_pp > -Inf
            mask[] <- FALSE
            mask[1] <- TRUE
            mask
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    p$initial_n[] <- 0
    p$initial_n_pp[] <- 0
    p$rr_pp[] <- 0.5
    p$cc_pp[] <- 0
    p$rr_pp[1] <- 0
    p <- setProcess(p, "resource",
                    prescribed_mask_fun = "stability_resource_input")
    result <- getStability(p, scope = "conditional", fixed = "consumers")
    expect_true(result$stable)
    expect_equal(result$n_active, length(p$w_full) - 1L)
    expect_equal(result$max_real_part, -0.5)
    expect_true(result$prescribed_masks$resource[1])
    expect_true(all(result$leading_eigenvectors$resource[1, ] == 0))
    expect_true(all(result$leading_eigenvectors$fish == 0))
    expect_equal(dim(result$leading_eigenvectors$resource),
                 c(length(p$w_full), 2L))
})

test_that("a complete scope with no dynamic coordinates is inapplicable", {
    p <- setResource(NS_params_small, resource_dynamics = "resource_constant",
                     balance = FALSE)
    for (getter in list(getStability, getDiscreteStability)) {
        result <- getter(p, scope = "conditional", fixed = "consumers")
        expect_identical(result$status, "not_applicable")
        expect_true(result$coverage$complete)
        expect_identical(result$n_active, 0L)
        expect_true(is.na(result$stable))
        expect_identical(result$prescribed, "resource")
        expect_equal(dim(result$leading_eigenvectors$resource),
                     c(length(p$w_full), 0L))
    }
})

test_that("boundary differencing is second order without negative probes", {
    seen <- numeric()
    ctx <- list(x0 = 0, scale = 1,
                check = function(...) NULL, describe = function(...) "trial")
    f <- function(x) {
        seen <<- c(seen, x)
        -2 * x + x^2
    }
    expect_equal(stability_jacobian(ctx, f, h = 0.2)[1, 1], -2)
    expect_equal(seen, c(0, 0.2, 0.4))
})

test_that("custom steps do not imply continuous equations or autonomy", {
    seen <- numeric()
    rlang::local_bindings(
        stability_custom_step = function(params, n_pp, t, dt, ...) {
            seen <<- c(seen, t)
            n_pp * exp(-dt)
        },
        stability_custom_assessment = function(params, n_pp) {
            data.frame(measure = "resource", value = sum(abs(n_pp)),
                       tolerance = 0, units = "density")
        },
        stability_custom_equation = function(params, n_pp, t) {
            seen <<- c(seen, t)
            -n_pp
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    p$initial_n_pp[] <- 0
    p <- setResource(p, resource_dynamics = "stability_custom_step",
                     balance = FALSE)
    p <- setProcess(p, "resource", role = "dynamic", autonomous = TRUE,
                    convergence_fun = "stability_custom_assessment")
    result <- getStability(p, scope = "conditional", fixed = "consumers", t = 7)
    expect_null(result$eigenvalues)
    expect_contains(result$coverage$unavailable$capability, "equation")
    discrete <- getDiscreteStability(p, scope = "conditional",
                                     fixed = "consumers", t = 7)
    expect_true(discrete$coverage$complete)
    expect_equal(discrete$spectral_radius, exp(-1), tolerance = 1e-8)
    expect_true(all(seen == 7))
    p <- setProcess(p, "resource", rate_fun = "stability_custom_equation")
    continuous <- getStability(p, scope = "conditional", fixed = "consumers",
                               t = 7)
    expect_equal(continuous$max_real_part, -1, tolerance = 1e-8)
    expect_identical(continuous$time, 7)
    expect_true(all(seen == 7))
    p <- setProcess(p, "resource", autonomous = FALSE)
    nonautonomous <- getStability(p, scope = "conditional", fixed = "consumers")
    expect_false(nonautonomous$autonomous)
    expect_null(nonautonomous$eigenvalues)
    expect_identical(nonautonomous$status, "inconclusive")
})

test_that("stability margins and provenance are explicit", {
    p <- NS_params_small
    p$initial_n[] <- 0
    p$initial_n_pp[] <- 0
    p$rr_pp[] <- 1e-9
    p$cc_pp[] <- 0
    marginal <- getStability(p, scope = "conditional", fixed = "consumers")
    strict <- getStability(p, scope = "conditional", fixed = "consumers",
                           stability_tol = 1e-10)
    expect_identical(marginal$status, "inconclusive")
    expect_true(strict$stable)
    expect_identical(strict$tolerance, 1e-10)
    expect_identical(attr(strict, "stability_key"),
                     stability_model_key(strict$params))
    strict$params$initial_n_pp[1] <- 1
    expect_false(identical(attr(strict, "stability_key"),
                           stability_model_key(strict$params)))
    expect_error(getStability(p, scope = "conditional"), "nonempty")
    expect_error(getStability(p, fixed = "resource"), "whole-model")
    expect_error(getStability(p, t = NA_real_), "finite")
    expect_error(getStability(p, stability_tol = -1), "stability_tol")
})

test_that("custom projection drivers cannot inherit the core Jacobian", {
    rlang::local_bindings(
        project.stability_driver = function(params, ...) NextMethod(),
        .env = .GlobalEnv
    )
    p <- recordExtension(NS_params_small, "stability_driver")
    class(p) <- c("stability_driver", class(p))
    for (owner in c("consumers", "resource")) {
        p <- setProcess(p, owner, role = "dynamic", autonomous = TRUE,
                        convergence_fun = NULL)
    }
    for (getter in list(getStability, getDiscreteStability)) {
        result <- getter(p)
        expect_identical(result$status, "inconclusive")
        expect_false(result$coverage$complete)
        expect_true(is.na(result$n_active))
        expect_true(any(grepl("projection driver",
                              result$coverage$unavailable$reason)))
    }
})

test_that("fixed contributions still require current temporal coverage", {
    rlang::local_bindings(
        stability_component_mort = function(params, ...) params$mu_b * 0,
        .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "forcing", initial_value = 1,
                      dynamics_fun = "constant_other",
                      mort_fun = "stability_component_mort")
    p <- setProcess(p, "consumers", autonomous = TRUE)
    p <- setProcess(p, "resource", autonomous = TRUE)
    missing <- getStability(p, scope = "conditional",
                            fixed = "component:forcing")
    expect_contains(missing$coverage$missing, "component:forcing")
    expect_null(missing$eigenvalues)
    excluded <- getStability(p, scope = "conditional",
                             fixed = c("consumers", "resource",
                                       "component:forcing"))
    expect_identical(excluded$status, "inconclusive")
    expect_true(is.na(excluded$n_active))
    p <- setComponent(p, "forcing", autonomous = TRUE)
    covered <- getStability(p, scope = "conditional",
                            fixed = "component:forcing")
    expect_true(covered$coverage$complete)
    expect_true(is.complex(covered$eigenvalues))
})

test_that("partial masks project the actual implicit Euler map", {
    rlang::local_bindings(
        stability_fish_input = function(params) {
            mask <- params$initial_n < 0
            j <- floor(mean(c(params$w_min_idx[1], support_top_idx(params)[1])))
            mask[1, j] <- TRUE
            mask
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    p$ext_diffusion[] <- outer(rep(1, nrow(p$initial_n)), 0.3 * p$w^2)
    p <- setProcess(p, "consumers",
                    prescribed_mask_fun = "stability_fish_input")
    ctx <- stability_context(p, scope = "conditional", fixed = "resource",
                             discrete = TRUE)
    dt <- 0.2
    mapped <- stability_step(ctx, dt)(ctx$x0)
    simulation <- project(p, t_max = dt, dt = dt, t_save = dt, method = "euler",
                          progress_bar = FALSE, info_level = 0)
    expect_equal(mapped, finalN(simulation)[ctx$active_idx])
    # Pinning excluded rows inside the implicit solve describes another map.
    rates <- getRates(p)
    coefs <- get_transport_coefs(
        p, n = p$initial_n, g = rates$e_growth, mu = rates$mort, dt = dt,
        recruitment_flux = rates$rdd, d = rates$diffusion,
        flux_limiter = flux_limiter_scheme(p)
    )
    mask <- ctx$prescribed_masks$consumers
    coefs$a[mask] <- 0
    coefs$b[mask] <- 1
    coefs$c[mask] <- 0
    coefs$S[mask] <- p$initial_n[mask]
    pinned <- project_n_loop(p$initial_n, coefs$a, coefs$b, coefs$c, coefs$S,
                             p$w_min_idx)
    expect_gt(max(abs(mapped - pinned[ctx$active_idx])), 1e-8)
    perturbed <- ctx$unpack(ctx$x0 * 1.01)
    expect_equal(perturbed$n[mask], p$initial_n[mask])
})

test_that("undefined logistic boundary equations produce partial coverage", {
    p <- NS_params_small
    p$initial_n[] <- 0
    p <- setResource(p, resource_dynamics = "resource_logistic",
                     balance = FALSE)
    p$initial_n_pp[] <- 0
    p$cc_pp[] <- 0
    p$rr_pp[] <- 1
    continuous <- getStability(p, scope = "conditional", fixed = "consumers")
    expect_identical(continuous$status, "inconclusive")
    expect_contains(continuous$coverage$missing, "resource")
    expect_null(continuous$eigenvalues)
    expect_true(any(grepl("capacity", continuous$coverage$unavailable$reason)))
    discrete <- getDiscreteStability(p, scope = "conditional",
                                     fixed = "consumers")
    expect_true(discrete$coverage$complete)
    expect_equal(discrete$n_active, length(p$w_full))
})

test_that("quadratic and mixed feedback keep genuine neutral coordinates", {
    rlang::local_bindings(
        stability_nonlinear_equation = function(params, n_pp) {
            derivative <- -n_pp
            derivative[1] <- -n_pp[1]^2
            derivative[2] <- n_pp[1] * n_pp[2]
            derivative
        },
        stability_nonlinear_step = function(params, n_pp, dt, ...) {
            n_pp + dt * stability_nonlinear_equation(params, n_pp)
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    p$initial_n_pp[] <- 0
    p <- setResource(p, resource_dynamics = "stability_nonlinear_step",
                     balance = FALSE)
    p <- setProcess(p, "resource", role = "dynamic", autonomous = TRUE,
                    rate_fun = "stability_nonlinear_equation",
                    convergence_fun = "builtin_resource_assessment")
    continuous <- getStability(p, scope = "conditional", fixed = "consumers")
    discrete <- getDiscreteStability(p, scope = "conditional",
                                     fixed = "consumers", dt = 0.1)
    expect_equal(continuous$n_active, length(p$w_full))
    expect_equal(discrete$n_active, length(p$w_full))
    expect_equal(Re(continuous$eigenvalues[1:2]), c(0, 0), tolerance = 1e-14)
    expect_equal(Re(discrete$discrete_eigenvalues[1:2]), c(1, 1))
    expect_identical(continuous$status, "inconclusive")
    expect_identical(discrete$status, "inconclusive")
})
