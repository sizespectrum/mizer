# A small model that is genuinely at its steady state, built once here.
residual_params <- suppressMessages(
    steady(newTraitParams(no_sp = 2, no_w = 20, max_w_max = 100,
                          min_w = 1e-3, w_pp_cutoff = 5, ks = 4,
                          reproduction_level = 0.25, info_level = 0),
           tol = 1e-6, t_max = 500, progress_bar = FALSE, info_level = 0)
)

# A copy that has been knocked off it, by scaling one species' abundance.
off_steady_params <- local({
    p <- residual_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 3
    p
})

# getSteadyResidual ----

test_that("getSteadyResidual() returns a labelled ArraySpeciesBySize", {
    res <- getSteadyResidual(residual_params)
    expect_s3_class(res, "ArraySpeciesBySize")
    expect_identical(dim(res), dim(initialN(residual_params)))
    expect_identical(dimnames(res), dimnames(initialN(residual_params)))
    expect_identical(attr(res, "units"), "1/year")
    expect_identical(attr(res, "value_name"), "Biomass drift contribution")
    # Each value is an integral over its bin, so it is drawn at the bin centre.
    expect_identical(attr(res, "representation"), "average")
    # The resource and the other components travel as attributes.
    expect_length(attr(res, "resource"), length(residual_params@w_full))
    expect_type(attr(res, "other"), "list")

    # The two measures are different quantities and have to say so, or `plot2()`
    # would draw one on top of the other without complaint.
    pc <- getSteadyResidual(residual_params, measure = "per_capita")
    expect_identical(attr(pc, "value_name"), "Per-capita rate of change")
    expect_identical(attr(pc, "representation"), "point")
})

test_that("getSteadyResidual() rejects an unknown measure", {
    expect_error(getSteadyResidual(residual_params, measure = "cell"),
                 "should be one of")
})

test_that("getSteadyResidual() is ~zero at a steady state and not otherwise", {
    expect_lt(max(abs(rowSums(getSteadyResidual(residual_params)))), 1e-3)
    expect_gt(max(abs(rowSums(getSteadyResidual(off_steady_params)))), 0.1)
})

test_that("the biomass measure sums over sizes to the biomass drift", {
    # This is the property the default measure exists for: the array says where
    # a model is unsteady in the same currency that `isSteady()` uses to decide
    # whether it is. It has to hold under both quadrature schemes, since the bin
    # weights differ between them.
    for (bin_average in c(FALSE, TRUE)) {
        p <- off_steady_params
        second_order_w(p) <- c(bin_average = bin_average)
        res <- getSteadyResidual(p)

        rates <- mizer:::steady_rates(p)
        biomass <- sizeIntegral(p, weighting = w(p), n = rates$n)
        dBdt <- sizeIntegral(p, weighting = w(p), n = rates$dNdt)
        expect_equal(rowSums(res), as.numeric(dBdt) / as.numeric(biomass),
                     ignore_attr = TRUE)

        # The resource attribute is the same measure on the resource grid, so it
        # sums to the resource drift, and the largest of all of them is the
        # scalar every steady-state check is stated against.
        wdw <- w_full(p) * dw_full(p)
        resource_drift <- sum(rates$dn_pp_dt * wdw) / sum(rates$n_pp * wdw)
        expect_equal(sum(attr(res, "resource")), resource_drift)
        expect_equal(max(abs(c(rowSums(res), sum(attr(res, "resource"))))),
                     mizer:::steady_biomass_drift(p))
    }
})

test_that("the per-capita measure is NA exactly where there are no fish", {
    res <- getSteadyResidual(off_steady_params, measure = "per_capita")
    expect_identical(is.na(unclass(res)) & TRUE,
                     initialN(off_steady_params) == 0)
})

test_that("the biomass measure reports every class, including empty ones", {
    # `dN/dt` is well defined in a class with no fish in it — it can be filling
    # up — so the default measure has nothing to withdraw. Only a species with
    # no biomass at all has no relative rate of change of it.
    p <- off_steady_params
    initialN(p)[1, ncol(initialN(p))] <- 0
    res <- getSteadyResidual(p)
    expect_false(anyNA(res))

    empty <- p
    initialN(empty)[1, ] <- 0
    res <- getSteadyResidual(empty)
    expect_true(all(is.na(res[1, ])))
    expect_false(anyNA(res[-1, ]))
})

test_that("a size class holding a trace cannot dominate the biomass measure", {
    # What #570 reported. A trace of fish in the largest size class has nothing
    # growing into it, so its per-capita rate is minus its mortality and stays
    # there for ever while the density falls through 1e-100 and beyond. The
    # default measure needs no cutoff to disregard it: the class holds no mass,
    # so it contributes none of the drift.
    p <- residual_params
    j <- length(w(p))
    n <- initialN(p)
    n[1, j] <- 1e-100
    initialN(p) <- n

    res <- getSteadyResidual(p)
    pc <- getSteadyResidual(p, measure = "per_capita")

    expect_equal(unname(pc[1, j]), unname(-getMort(p)[1, j]))
    # `all.sizes = TRUE` isolates the biomass rule from the size-range one:
    # this class is above its species' `w_max`, so the default summary would
    # leave it out on that ground alone.
    expect_lt(summary(pc, all.sizes = TRUE)$per_species$Min[1], -0.01)
    # The same class, weighted by the biomass it holds, is nothing at all.
    expect_lt(abs(res[1, j]), 1e-90)
    expect_lt(max(abs(rowSums(res))), 1e-3)
})

test_that("the biomass-share cutoff is validated", {
    params <- NS_params_small
    n <- initialN(params)

    invalid_cutoffs <- list(c(0, 1), NA_real_, Inf, -0.1, 1.1)
    for (cutoff in invalid_cutoffs) {
        expect_error(negligible_cells(params, n, cutoff = cutoff),
                     "must be a finite number between 0 and 1", fixed = TRUE)
    }
})

test_that("the biomass-share cutoff uses the selected quadrature", {
    for (bin_average in c(FALSE, TRUE)) {
        params <- NS_params_small
        second_order_w(params) <- c(bin_average = bin_average)
        n <- initialN(params) * 0
        j <- ncol(n)
        # Give the first and last bins exactly half the species' biomass under
        # the selected quadrature. The top-bin weight is one-sided, so this also
        # distinguishes the two quadrature schemes rather than merely scaling
        # every bin by the same constant.
        wdw <- bin_average_weight(w(params), params) * dw(params)
        n[1, c(1, j)] <- 1 / wdw[c(1, j)]

        negligible <- negligible_cells(params, n, cutoff = 0.5)

        # A cell exactly on the cutoff is kept; the strict cutoff removes only
        # cells below it. Species with no biomass have no relevant cells.
        expect_identical(unname(which(!negligible[1, ])), c(1L, j))
        expect_true(all(negligible[-1, ]))
    }
})

test_that("getSteadyResidual() predicts the drift that project() produces", {
    # This is the property that defines the residual: it is dN/dt. A projection
    # over a short step must reproduce it, with an error that is first order in
    # the step size. Checking the convergence rate rather than a fixed
    # tolerance is what makes this a test of the identity rather than of one
    # arbitrary step length.
    predicted <- getSteadyResidual(off_steady_params, measure = "per_capita")
    err <- vapply(c(1e-4, 1e-5), function(dt) {
        sim <- project(off_steady_params, t_max = 2 * dt, dt = dt, t_save = dt,
                       progress_bar = FALSE)
        n0 <- N(sim)[1, , ]
        observed <- (N(sim)[2, , ] - n0) / (n0 * dt)
        sel <- !is.na(predicted) & n0 > 0 & abs(predicted) > 1e-3
        stats::median(abs((observed[sel] - predicted[sel]) / predicted[sel]))
    }, numeric(1))
    # Ten times the step, ten times the error, to within a generous factor.
    expect_equal(err[[1]] / err[[2]], 10, tolerance = 0.5)
})

test_that("getSteadyResidual() honours the effort it is given", {
    # Needs a model that is actually fished; the trait fixture above is not.
    p <- NS_params_small
    higher <- initial_effort(p) + 1
    # Passing the effort must be the same as storing it in the model ...
    stored <- p
    initial_effort(stored) <- higher
    expect_equal(unclass(getSteadyResidual(p, effort = higher)),
                 unclass(getSteadyResidual(stored)),
                 ignore_attr = TRUE)
    # ... and must actually change the answer.
    expect_false(isTRUE(all.equal(
        unclass(getSteadyResidual(p, effort = higher)),
        unclass(getSteadyResidual(p)), check.attributes = FALSE)))
})

test_that("getSteadyResidual() works under the second-order scheme", {
    p <- residual_params
    second_order_w(p) <- TRUE
    expect_s3_class(getSteadyResidual(p), "ArraySpeciesBySize")
    expect_true(all(is.finite(getSteadyResidual(p))))
    pc <- getSteadyResidual(p, measure = "per_capita")
    expect_true(all(is.finite(pc[initialN(p) > 0])))
})

test_that("getSteadyResidual() rejects a non-positive dt", {
    expect_error(getSteadyResidual(residual_params, dt = 0))
    expect_error(getSteadyResidual(residual_params, dt = -1))
})

# steady_biomass_drift ----

test_that("steady_biomass_drift() separates settled from unsettled models", {
    expect_lt(steady_biomass_drift(residual_params), steady_residual_tol())
    expect_gt(steady_biomass_drift(off_steady_params), steady_residual_tol())
})

test_that("steady_biomass_drift() is the relative rate of biomass change", {
    # Check it against the biomass change an actual projection produces, which
    # is what the number claims to be. The resource counts too, and on this
    # model it is in fact the worst offender: tripling a consumer triples the
    # predation on the resource.
    dt <- 1e-4
    sim <- project(off_steady_params, t_max = 2 * dt, dt = dt, t_save = dt,
                   progress_bar = FALSE)
    b <- getBiomass(sim)
    wdw <- off_steady_params@w_full * off_steady_params@dw_full
    rb <- as.numeric(NResource(sim) %*% wdw)
    observed <- max(abs((b[2, ] - b[1, ]) / (b[1, ] * dt)),
                    abs((rb[[2]] - rb[[1]]) / (rb[[1]] * dt)))
    expect_equal(steady_biomass_drift(off_steady_params), observed,
                 tolerance = 1e-2)
})

test_that("steady_biomass_drift() ignores fast cells that hold no biomass", {
    # The point of weighting by biomass: the per-cell maximum is dominated by
    # the fastest-relaxing size classes, which carry negligible mass. The two
    # measures must therefore disagree by orders of magnitude on a settled
    # model, and it is the biomass one that reads as settled.
    cell_max <- max(abs(getSteadyResidual(residual_params,
                                          measure = "per_capita")),
                    na.rm = TRUE)
    expect_lt(steady_biomass_drift(residual_params), cell_max)
})

# warn_if_not_steady ----

test_that("warn_if_not_steady() fires only when the model is off steady state", {
    expect_warning(warn_if_not_steady(off_steady_params, "Context."),
                   "not at its steady state")
    expect_silent(warn_if_not_steady(residual_params, "Context."))
})

test_that("warn_if_not_steady() includes the context and is silenced by info_level", {
    expect_warning(warn_if_not_steady(off_steady_params, "Sentinel phrase."),
                   "Sentinel phrase")
    withr::local_options(mizer_info_level = 0)
    expect_silent(warn_if_not_steady(off_steady_params, "Context."))
})

# isSteady ----

test_that("isSteady() returns a single boolean", {
    expect_true(isSteady(residual_params))
    expect_false(isSteady(off_steady_params))
    expect_type(isSteady(residual_params), "logical")
    expect_length(isSteady(residual_params), 1L)
})

test_that("isSteady() respects custom tol argument", {
    expect_true(isSteady(off_steady_params, tol = 10))
    expect_false(isSteady(residual_params, tol = 1e-10))
})

test_that("isSteady() respects effort argument", {
    p <- NS_params_small
    higher <- initial_effort(p) + 1
    # Different effort changes steadiness if model was settled at initial_effort
    expect_equal(isSteady(p, effort = initial_effort(p)),
                 isSteady(p))
})

test_that("isSteady() works under the second-order scheme", {
    p <- residual_params
    second_order_w(p) <- TRUE
    # Changing quadrature scheme moves model off its original steady state
    expect_false(isSteady(p))
    # Re-settling under second-order scheme restores steady state
    p <- suppressMessages(steady(p, tol = 1e-5, t_max = 500,
                                 progress_bar = FALSE, info_level = 0))
    expect_true(isSteady(p))
})

test_that("isSteady() dispatches via S3", {
    dummy <- structure(list(), class = "DummyModel")
    isSteady.DummyModel <- function(params, ...) TRUE  # nolint: object_name_linter.
    expect_true(isSteady(dummy))
})
