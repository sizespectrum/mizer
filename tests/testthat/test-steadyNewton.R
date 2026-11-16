# Tests for the internals of the direct (Newton) steady-state solver.
# Its exported entry points are tested in test-steadyState.R; stability
# analysis is tested in test-stability.R.

# A steadied model is a good starting guess for the root finder.
delayedAssign("p_steady",
               tuneSteadyState(NS_params_small, t_max = 50,
                               progress_bar = FALSE))

test_that("Newton residuals reject non-finite trials before evaluating rates", {
    p <- NS_params_small
    second_order_w(p) <- TRUE
    for (resource in c("fixed", "solve")) {
        active <- steady_active_set(p, resource = resource)
        n0 <- positive_initial_guess(p@initial_n, p@w_min_idx, active$w_top)
        f <- steady_state_residual(p, getRDD(p), p@initial_n_other,
                                   p@initial_effort, active, n0 = n0)
        x <- log(n0[active$mask])
        if (resource == "solve") {
            x <- c(x, log(p@initial_n_pp[active$mask_pp]))
        }
        expect_true(all(is.finite(f(x))))
        for (bad in c(NA_real_, Inf, -Inf, 1000, -1000)) {
            for (index in c(1L, length(x))) {
                trial <- x
                trial[index] <- bad
                expect_equal(f(trial), rep(Inf, length(x)))
            }
        }
    }
})

test_that("Newton safeguards do not swallow custom rate-function errors", {
    p <- NS_params_small
    rdd <- getRDD(p)
    assign("newton_test_growth_error", function(...) {
        stop("custom growth failed")
    }, envir = .GlobalEnv)
    withr::defer(rm("newton_test_growth_error", envir = .GlobalEnv))
    # Register directly to exercise a runtime failure rather than the setter's
    # own test call of the custom function.
    p@rates_funcs$EGrowth <- "newton_test_growth_error"
    active <- steady_active_set(p, resource = "fixed")
    n0 <- positive_initial_guess(p@initial_n, p@w_min_idx, active$w_top)
    f <- steady_state_residual(p, rdd, p@initial_n_other,
                               p@initial_effort, active, n0 = n0)
    expect_equal(f(rep(1000, active$n_fish_active)),
                 rep(Inf, active$n_fish_active))
    expect_error(f(log(n0[active$mask])), "custom growth failed")
    expect_error(newton_frozen_guess(p, active, n0, p@initial_effort),
                 "custom growth failed")
})

test_that("the frozen Newton guess keeps diffusion under custom growth", {
    p <- NS_params_small
    second_order_w(p) <- TRUE
    assign("newton_test_zero_growth", function(params, ...) {
        params@initial_n * 0
    }, envir = .GlobalEnv)
    withr::defer(rm("newton_test_zero_growth", envir = .GlobalEnv))
    p <- setRateFunction(p, "EGrowth", "newton_test_zero_growth")
    ext_diffusion(p) <- 0.1 * outer(p@species_params$n, p@w,
                                    function(n, w) w^(n + 1))
    active <- steady_active_set(p, resource = "fixed")
    n0 <- positive_initial_guess(p@initial_n, p@w_min_idx, active$w_top)
    seed <- newton_frozen_guess(p, active, n0, p@initial_effort)
    expect_true(all(seed[active$mask] > 0))
    expect_true(all(seed[!active$mask] == 0))
    expect_equal(active$w_top, unname(support_top_idx(p)))
    # Diffusion feeds cells even though the supplied growth is exactly zero.
    next_cell <- cbind(seq_len(nrow(seed)), p@w_min_idx + 1L)
    expect_true(all(seed[next_cell] > apply(n0, 1, max) * 1e-15))
})

test_that("bounded Newton steps solve the North Sea zero-growth tail", {
    skip_unless_experimental()
    # The full data are needed for the reported Cod cell at the final grid
    # boundary; the smaller fixture does not reproduce its near-singular column.
    p <- NS_params
    second_order_w(p) <- TRUE
    top <- support_top_idx(p)
    sol <- newton_steady_state(p, p@initial_effort, getRDD(p),
                               resource = "fixed")
    expect_lt(sol$residual, 1e-6)
    expect_true(all(is.finite(sol$n)))
    cod <- match("Cod", p@species_params$species)
    expect_equal(unname(sol$n[cod, top[cod]]), 0)
    expect_equal(support_top_idx(p), top)
})

test_that("the Newton retry solves a model whose w_max exceeds w_inf", {
    skip_unless_experimental()
    sp <- data.frame(species = "test", w_inf = 1000, w_mat = 100)
    p <- newMultispeciesParams(sp, no_w = 50, second_order_w = TRUE,
                               info_level = 0)
    expect_equal(unname(p@species_params$w_max), 1500)
    top <- support_top_idx(p)
    sol <- newton_steady_state(p, p@initial_effort, getRDD(p),
                               resource = "fixed")
    expect_lt(sol$residual, 1e-6)
    expect_true(all(is.finite(sol$n)))
    expect_equal(support_top_idx(p), top)
})

test_that("the Newton retry retains a diffusion-fed North Sea tail", {
    skip_unless_experimental()
    p <- NS_params
    second_order_w(p) <- TRUE
    ext_diffusion(p) <- 0.1 * outer(p@species_params$n, p@w,
                                    function(n, w) w^(n + 1))
    top <- support_top_idx(p)
    sol <- newton_steady_state(p, p@initial_effort, getRDD(p),
                               resource = "fixed")
    expect_lt(sol$residual, 1e-6)
    cod <- match("Cod", p@species_params$species)
    expect_gt(sol$n[cod, top[cod]], max(p@initial_n[cod, ]) * 1e-15)
    expect_true(all(is.finite(sol$n)))
    expect_equal(support_top_idx(p), top)
})

test_that("Newton retries unmet residuals and retains the better result", {
    skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")),
                "requires testthat::local_mocked_bindings()")
    p <- NS_params_small
    # Exercise the root finder's two distinct stopping criteria: a small step
    # (termcd 2) does not establish that the residual tolerance has been met.
    for (errors in list(c(1e-8), c(1, 0.1), c(0.1, 1))) {
        calls <- 0L
        seeds <- 0L
        budgets <- numeric(0)
        solver_stub <- function(x, fn, method, global, control) {
            calls <<- calls + 1L
            budgets <<- c(budgets, control$maxit)
            list(x = x, fvec = rep(errors[calls], length(x)),
                 iter = 3L, termcd = 2L, message = "small step")
        }
        local_mocked_bindings(nleqslv = solver_stub, .package = "nleqslv")
        guess_stub <- function(params, active, n0, effort, rates_fns) {
            seeds <<- seeds + 1L
            n0
        }
        local_mocked_bindings(newton_frozen_guess = guess_stub)
        sol <- newton_steady_state(p, p@initial_effort, getRDD(p),
                                   resource = "fixed", maxit = 10)
        expect_equal(calls, length(errors))
        expect_equal(seeds, length(errors) - 1L)
        expect_equal(sol$residual, min(errors))
        expect_equal(budgets, if (length(errors) == 1L) 10 else c(10, 7))
    }
})

test_that("a failed frozen-rate retry keeps the original Newton result", {
    skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")),
                "requires testthat::local_mocked_bindings()")
    p <- NS_params_small
    solver_stub <- function(x, fn, method, global, control) {
        list(x = x, fvec = rep(0.1, length(x)), iter = 1L,
             termcd = 5L, message = "ill-conditioned Jacobian")
    }
    local_mocked_bindings(nleqslv = solver_stub, .package = "nleqslv")
    local_mocked_bindings(newton_frozen_guess = function(...) {
        stop(structure(list(message = "frozen solve did not converge",
                            call = NULL),
                       class = c("mizer_steady_nonconvergence", "error",
                                 "condition")))
    })
    sol <- newton_steady_state(p, p@initial_effort, getRDD(p),
                               resource = "fixed")
    expect_equal(sol$residual, 0.1)
    expect_equal(sol$termcd, 5L)
})

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

test_that("consumer residuals use final rates from complete providers", {
    assign("newton_complete_rates", function(params, ...) {
        rates <- mizerRates(params, ...)
        rates$e_growth <- 0.8 * rates$e_growth
        rates$mort <- rates$mort + 0.2
        rates$diffusion <- rates$diffusion + 0.01
        rates$rdd <- 1.7 * rates$rdd
        rates
    }, envir = .GlobalEnv)
    withr::defer(rm("newton_complete_rates", envir = .GlobalEnv))
    p <- NS_params_small
    p$rates_funcs$Rates <- "newton_complete_rates"
    for (scheme in list(FALSE, TRUE)) {
        second_order_w(p) <- scheme
        rates <- getRates(p)
        expected <- consumer_residual(p, n = p$initial_n,
                                      n_pp = p$initial_n_pp,
                                      n_other = p$initial_n_other,
                                      effort = p$initial_effort,
                                      rates = rates, rdd = rates$rdd)
        actual <- consumer_residual(p, n = p$initial_n,
                                    n_pp = p$initial_n_pp,
                                    n_other = p$initial_n_other,
                                    effort = p$initial_effort)
        expect_equal(actual, expected, ignore_attr = TRUE)
        partial <- consumer_residual(
            p, n = p$initial_n, n_pp = p$initial_n_pp,
            n_other = p$initial_n_other, effort = p$initial_effort,
            rates = rates[c("e_growth", "mort", "diffusion")]
        )
        expect_equal(partial, expected, ignore_attr = TRUE)
    }
})

test_that("Newton shares complete rates across its equations", {
    p <- NS_params_small
    active <- steady_active_set(p, resource = "solve")
    n0 <- positive_initial_guess(p$initial_n, p$w_min_idx, active$w_top)
    target <- getRates(p)$rdd
    functions <- projectRateFunctions(p)
    baseline <- steady_state_residual(p, target, p$initial_n_other,
                                      p$initial_effort, active, n0 = n0,
                                      rates_fns = functions)
    original <- functions$Rates
    observation <- new.env(parent = emptyenv())
    observation$calls <- 0L
    functions$Rates <- function(...) {
        observation$calls <- observation$calls + 1L
        rates <- original(...)
        rates$mort <- rates$mort + 0.2
        rates$resource_mort <- rates$resource_mort + 0.4
        rates$rdd <- 1000 * rates$rdd
        rates
    }
    changed <- steady_state_residual(p, target, p$initial_n_other,
                                     p$initial_effort, active, n0 = n0,
                                     rates_fns = functions)
    x <- log(c(n0[active$mask], p$initial_n_pp[active$mask_pp]))
    difference <- changed(x) - baseline(x)
    expect_equal(observation$calls, 1L)
    expect_equal(unname(difference[seq_len(active$n_fish_active)]),
                 rep(0.2, active$n_fish_active), tolerance = 1e-7)
    expect_equal(unname(tail(difference, sum(active$mask_pp))),
                 rep(-0.4, sum(active$mask_pp)), tolerance = 1e-7)
})

test_that("Newton refuses consumer dynamics without its equation operator", {
    method <- "project_simple.operatorTest"
    assign(method, function(params, ...) {
        NextMethod()
    }, envir = .GlobalEnv)
    assign("operator_assessment", function(...) {
        data.frame(measure = "residual", value = 0, tolerance = 1,
                   units = "1/year")
    }, envir = .GlobalEnv)
    withr::defer(rm("project_simple.operatorTest", "operator_assessment",
                    envir = .GlobalEnv))
    p <- recordExtension(NS_params_small, "operatorTest")
    class(p) <- c("operatorTest", class(p))
    for (owner in c("consumers", "resource")) {
        p <- setProcess(p, owner, role = "dynamic",
                        convergence_fun = "operator_assessment")
    }
    expect_false(builtin_consumer_projection(p))
    expect_true(builtin_consumer_projection(NS_params_small))
    for (finder in list(findSteadyState, tuneSteadyState)) {
        expect_error(finder(p, solver = "newton", info_level = 0),
                     "built-in consumer transport equations")
    }
})

test_that("Newton retains a coupled root with state-dependent resource rates", {
    skip_unless_experimental()
    p <- p_steady
    reference_rates <- getRates(p)
    feedback <- function(params, n_pp, ...) {
        rates <- reference_rates
        scale <- ifelse(params$cc_pp > 0, params$cc_pp, 1)
        rates$resource_mort <- 100 * params$rr_pp * (n_pp / scale)^2
        rates
    }
    assign("newton_resource_feedback", feedback, envir = .GlobalEnv)
    withr::defer(rm("newton_resource_feedback", envir = .GlobalEnv))
    p$rates_funcs$Rates <- "newton_resource_feedback"
    sol <- newton_steady_state(p, effort = p$initial_effort,
                               rdd_const = reference_rates$rdd,
                               resource = "solve")
    active <- p$cc_pp > 0 & p$rr_pp > 0
    ratio <- sol$n_pp[active] / p$cc_pp[active]
    # The resource root satisfies x + 100*x^3 = 1. A fixed-point iteration
    # x <- 1/(1 + 100*x^2) is unstable there and must not replace this root.
    expect_equal(ratio + 100 * ratio^3, rep(1, sum(active)),
                 tolerance = 1e-6, ignore_attr = TRUE)
    p$initial_n[] <- sol$n
    p$initial_n_pp[] <- sol$n_pp
    rates <- getRates(p)
    residual <- p$rr_pp * (p$cc_pp - p$initial_n_pp) -
        rates$resource_mort * p$initial_n_pp
    scale <- p$rr_pp[active] * p$cc_pp[active]
    expect_lt(max(abs(residual[active] / scale)), 1e-6)
    expect_true(getSteadyAssessment(p)$verdict)
})
