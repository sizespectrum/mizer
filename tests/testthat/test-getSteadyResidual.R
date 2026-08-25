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
    expect_identical(attr(res, "value_name"), "Steady-state residual")
    # The resource and the other components travel as attributes.
    expect_length(attr(res, "resource"), length(residual_params@w_full))
    expect_type(attr(res, "other"), "list")
})

test_that("getSteadyResidual() is ~zero at a steady state and not otherwise", {
    expect_lt(max(abs(getSteadyResidual(residual_params)), na.rm = TRUE), 1e-3)
    expect_gt(max(abs(getSteadyResidual(off_steady_params)), na.rm = TRUE), 0.1)
})

test_that("getSteadyResidual() is NA exactly where there are no fish", {
    res <- getSteadyResidual(off_steady_params)
    expect_identical(is.na(unclass(res)) & TRUE,
                     initialN(off_steady_params) == 0)
})

test_that("getSteadyResidual() reports nothing for a negligible size class", {
    # A settled model, where every genuine cell rate is below 1e-3, so that the
    # trace added below is the only thing that could dominate a summary.
    p <- residual_params
    j <- length(w(p))
    n <- initialN(p)
    # A trace of fish in the largest size class. Nothing grows into it, so its
    # per-capita rate is minus its mortality and stays there for ever while the
    # density falls through 1e-100 and beyond. Reporting it would let a class
    # holding no mass dominate the summary and the plot.
    n[1, j] <- 1e-100
    initialN(p) <- n

    full <- getSteadyResidual(p, biomass_share_cutoff = 0)
    trimmed <- getSteadyResidual(p)

    expect_false(is.na(full[1, j]))
    expect_equal(unname(full[1, j]), unname(-getMort(p)[1, j]))
    expect_true(is.na(trimmed[1, j]))
    # Nothing else moved: exactly the one class was withdrawn.
    expect_identical(sum(is.na(unclass(trimmed))),
                     sum(is.na(unclass(full))) + 1L)
    # This is what #570 reported: the trace dominates the default summary of a
    # model that is otherwise settled to within 1e-3.
    expect_lt(summary(full)$per_species$Min[1], -0.01)
    expect_gt(summary(trimmed)$per_species$Min[1], -1e-3)
})

test_that("getSteadyResidual() is unchanged where every class holds fish", {
    # The cutoff only ever removes classes that hold no mass, so on a model
    # without such a class it changes nothing at all.
    expect_identical(unclass(getSteadyResidual(off_steady_params)),
                     unclass(getSteadyResidual(off_steady_params,
                                               biomass_share_cutoff = 0)))
})

test_that("getSteadyResidual() predicts the drift that project() produces", {
    # This is the property that defines the residual: it is dN/dt. A projection
    # over a short step must reproduce it, with an error that is first order in
    # the step size. Checking the convergence rate rather than a fixed
    # tolerance is what makes this a test of the identity rather than of one
    # arbitrary step length.
    predicted <- getSteadyResidual(off_steady_params)
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
    res <- getSteadyResidual(p)
    expect_s3_class(res, "ArraySpeciesBySize")
    expect_true(all(is.finite(res[initialN(p) > 0])))
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
    cell_max <- max(abs(getSteadyResidual(residual_params)), na.rm = TRUE)
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
