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

# Component dynamics for the tests of what the drift does and does not cover.
# The state is grown by a fixed exponential factor each step, so its relative
# rate of change is a known constant, independent of `dt` to first order.
e <- globalenv()
e$drifting_component <- function(params, n_other, component, dt, ...) {
    n_other[[component]] * exp(2 * dt)
}
e$settled_component <- function(params, n_other, component, ...) {
    n_other[[component]]
}
e$opaque_component <- function(params, n_other, component, ...) {
    n_other[[component]]
}
e$drifting_component_assessment <- function(...) {
    data.frame(measure = "relative_change", value = 2,
               tolerance = 0.05, units = "1/year")
}
e$settled_component_assessment <- function(...) {
    data.frame(measure = "relative_change", value = 0,
               tolerance = 0.05, units = "1/year")
}

# A settled model carrying one component that never settles.
component_params <- setComponent(
    residual_params, "drifter", initial_value = c(1, 2, 3),
    dynamics_fun = "drifting_component",
    convergence_fun = "drifting_component_assessment"
)

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

test_that("empty-tail drift agrees with short projected mass changes", {
    for (flux in c("upwind", "van_leer")) {
        for (bin_average in c(FALSE, TRUE)) {
            p <- NS_params_small
            second_order_w(p) <- c(flux = flux, bin_average = bin_average)
            ext_diffusion(p) <- 0.1 * outer(p@species_params$n, p@w,
                                            function(n, w) w^(n + 1))
            outside <- col(p@initial_n) > support_top_idx(p)
            p@initial_n[outside] <- 0
            biomass <- as.numeric(getBiomass(p))
            numbers <- as.numeric(getN(p))
            predicted_b <- rowSums(getSteadyResidual(p))
            predicted_n <- drop(steady_rates(p)$dNdt %*% p@dw) / numbers
            error <- vapply(c(1e-6, 1e-7), function(dt) {
                sim <- project(p, t_max = dt, dt = dt, t_save = dt,
                               method = "euler", progress_bar = FALSE,
                               info_level = 0)
                next_p <- finalParams(sim)
                observed_b <- (as.numeric(getBiomass(next_p)) - biomass) /
                    (dt * biomass)
                observed_n <- (as.numeric(getN(next_p)) - numbers) /
                    (dt * numbers)
                expect_true(all(initialN(next_p)[outside] == 0))
                c(biomass = max(abs(observed_b - predicted_b)),
                  numbers = max(abs(observed_n - predicted_n)))
            }, numeric(2))
            # Time-step error shrinks; spurious tail influx instead leaves a
            # nonzero discrepancy as dt tends to zero.
            expect_true(all(error[, 1] / error[, 2] > 8))
            expect_true(all(error[, 1] / error[, 2] < 12))
            expect_lt(error["biomass", 2], 2e-5)
            expect_lt(error["numbers", 2], 1.5e-4)
        }
    }
})

test_that("tail correction preserves occupied and retained equations", {
    for (flux in c("upwind", "van_leer", "centred")) {
        p <- NS_params_small
        second_order_w(p) <- c(flux = flux, bin_average = TRUE)
        ext_diffusion(p) <- 0.1 * outer(p@species_params$n, p@w,
                                        function(n, w) w^(n + 1))
        outside <- col(p@initial_n) > support_top_idx(p)
        p@initial_n[outside] <- 0
        inside_empty <- support_top_idx(p)[1] - 1L
        p@initial_n[1, inside_empty] <- 0
        p@initial_n[1, ncol(p@initial_n)] <- 1e-100
        n <- p@initial_n
        raw <- -consumer_residual(p, n, p@initial_n_pp, p@initial_n_other,
                                  p@initial_effort)
        realised <- steady_rates(p)$dNdt
        empty_tail <- outside & n == 0
        expect_true(any(raw[empty_tail] > 0))
        expect_true(all(realised[empty_tail] == 0))
        expect_equal(realised[!empty_tail], raw[!empty_tail])
        expect_gt(realised[1, inside_empty], 0)
        expect_equal(realised[1, ncol(n)], raw[1, ncol(n)])
        # Existing populated tails, however small, keep their per-capita rate.
        pc <- getSteadyResidual(p, measure = "per_capita")
        expect_equal(unname(pc[1, ncol(n)]),
                     unname(raw[1, ncol(n)] / n[1, ncol(n)]))
        expect_true(all(is.na(pc[empty_tail])))
    }
})

test_that("tail correction uses global diffusion and keeps legacy flow", {
    assign("steady_tail_test_growth", function(params, ...) {
        params@other_params$tail_growth
    }, envir = .GlobalEnv)
    assign("steady_tail_test_diffusion", function(params, ...) {
        params@other_params$tail_diffusion
    }, envir = .GlobalEnv)
    withr::defer(rm(list = c("steady_tail_test_growth",
                             "steady_tail_test_diffusion"), envir = .GlobalEnv))
    p <- NS_params_small
    second_order_w(p) <- FALSE
    top <- support_top_idx(p)[1]
    p@initial_n[] <- 0
    p@initial_n[1, top] <- 1
    p@other_params$tail_growth <- p@initial_n
    p@other_params$tail_diffusion <- p@initial_n * 0
    p <- setRateFunction(p, "EGrowth", "steady_tail_test_growth")
    p <- setRateFunction(p, "Diffusion", "steady_tail_test_diffusion")
    expect_true(all(p@ext_diffusion == 0))

    # With no diffusion, the legacy Euler step actually fills this empty tail.
    no_diffusion <- steady_rates(p)$dNdt
    expect_gt(no_diffusion[1, top + 1L], 0)
    dt <- 1e-6
    sim <- project(p, t_max = dt, dt = dt, t_save = dt, method = "euler",
                   progress_bar = FALSE, info_level = 0)
    expect_equal(finalN(sim)[1, top + 1L] / dt,
                 no_diffusion[1, top + 1L],
                 tolerance = 1e-5, ignore_attr = TRUE)

    # Diffusion on Herring selects the hard-boundary solver for Sprat too.
    # It comes entirely from the custom rate, not from a parameter slot.
    p@other_params$tail_diffusion[2, top] <- 1
    global_diffusion <- steady_rates(p)$dNdt
    expect_equal(unname(global_diffusion[1, top + 1L]), 0)
    expect_equal(global_diffusion[1, top], no_diffusion[1, top])
    sim <- project(p, t_max = dt, dt = dt, t_save = dt, method = "euler",
                   progress_bar = FALSE, info_level = 0)
    expect_equal(finalN(sim)[1, top + 1L], 0)

    # A non-upwind scheme enforces the same boundary even without diffusion.
    for (flux in c("van_leer", "centred")) {
        second_order_w(p) <- c(flux = flux)
        top <- support_top_idx(p)[1]
        p@initial_n[] <- 0
        p@initial_n[1, top] <- 1
        p@other_params$tail_growth[] <- 1
        p@other_params$tail_diffusion[] <- 0
        raw <- -consumer_residual(p, p@initial_n, p@initial_n_pp,
                                  p@initial_n_other, p@initial_effort)
        expect_gt(raw[1, top + 1L], 0)
        expect_equal(unname(steady_rates(p)$dNdt[1, top + 1L]), 0)
    }
})

test_that("tail correction keeps recruitment and outgoing number losses", {
    assign("steady_tail_test_zero", function(params, ...) {
        params@initial_n * 0
    }, envir = .GlobalEnv)
    assign("steady_tail_test_births", function(params, ...) {
        rep(0.25, nrow(params@initial_n))
    }, envir = .GlobalEnv)
    withr::defer(rm(list = c("steady_tail_test_zero",
                             "steady_tail_test_births"), envir = .GlobalEnv))
    for (flux in c("upwind", "van_leer")) {
        for (bin_average in c(FALSE, TRUE)) {
            p <- NS_params_small
            second_order_w(p) <- c(flux = flux, bin_average = bin_average)
            top <- support_top_idx(p)
            cells <- cbind(seq_len(nrow(p@initial_n)), top)
            p@initial_n[] <- 0
            p@initial_n[cells] <- 1 / p@dw[top]
            ext_diffusion(p) <- outer(rep(0.1, nrow(p@initial_n)), p@w^2)
            p <- setRateFunction(p, "EGrowth", "steady_tail_test_zero")
            p <- setRateFunction(p, "Mort", "steady_tail_test_zero")
            p <- setRateFunction(p, "RDD", "steady_tail_test_births")
            diffusion <- getDiffusion(p)
            # Integrating the PDE leaves births minus the outgoing upper-face
            # flux. This includes Cod's native terminal bin, which has no next
            # bin on the grid and must not be removed from the diagnostic.
            face_width <- if (flux == "upwind") {
                p@dw[top]
            } else {
                log(p@w[2] / p@w[1]) * p@w[top] * p@w[2] / p@w[1]
            }
            outflow <- diffusion[cells] * p@initial_n[cells] / (2 * face_width)
            rates <- steady_rates(p)$dNdt
            expect_equal(drop(rates %*% p@dw), 0.25 - outflow,
                         ignore_attr = TRUE)
            expect_true(all(outflow > 0))
            expect_equal(top[3], ncol(p@initial_n), ignore_attr = TRUE)
            eggs <- cbind(seq_len(nrow(rates)), p@w_min_idx)
            expect_true(all(rates[eggs] > 0))
            dt <- 1e-5
            sim <- project(p, t_max = dt, dt = dt, t_save = dt,
                           method = "euler", progress_bar = FALSE,
                           info_level = 0)
            observed <- (as.numeric(getN(finalParams(sim))) -
                             as.numeric(getN(p))) / dt
            expect_equal(observed, 0.25 - outflow,
                         tolerance = 1e-5, ignore_attr = TRUE)
        }
    }
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
        unclass(getSteadyResidual(p)), check.attributes = FALSE
    )))
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

# steady_drift_report ----

test_that("steady_drift_report() names the state variable responsible", {
    # Reporting the number without the name is what sent the reporter of #589
    # to the wrong place. Whichever variable wins, the name has to be the one
    # that achieves the reported drift.
    for (p in list(residual_params, off_steady_params)) {
        report <- steady_drift_report(p)
        expect_identical(report$drift, steady_biomass_drift(p))
        species <- abs(rowSums(getSteadyResidual(p)))
        wdw <- w_full(p) * dw_full(p)
        rates <- steady_rates(p)
        resource <- abs(sum(rates$dn_pp_dt * wdw) / sum(rates$n_pp * wdw))
        if (report$is_resource) {
            expect_equal(report$drift, resource)
            expect_gte(resource, max(species))
        } else {
            expect_equal(report$drift, species[[report$variable]])
            expect_gte(species[[report$variable]], resource)
        }
    }
})

test_that("steady_drift_report() names a species when a species is worst", {
    # Pinning the resource takes it out of the running, so a species must win.
    # Tripling one consumer also triples the predation on the resource, which
    # is why the resource is otherwise the worst offender on this model.
    p <- off_steady_params
    resource_dynamics(p) <- "resource_constant"
    report <- steady_drift_report(p)
    expect_false(report$is_resource)
    expect_identical(report$variable,
                     as.character(species_params(p)$species[[1]]))
    expect_identical(steady_variable_txt(report),
                     paste0("`", report$variable, "`"))
    expect_identical(steady_variable_txt(report, quote = FALSE),
                     report$variable)
})

test_that("steady_drift_report() reports the resource as the resource", {
    p <- residual_params
    initialNResource(p) <- initialNResource(p) * 3
    report <- steady_drift_report(p)
    expect_true(report$is_resource)
    expect_identical(steady_variable_txt(report), "the resource")
    expect_identical(steady_variable_txt(report, quote = FALSE), "the resource")
})

test_that("steady_drift_report() covers only components that can move", {
    # A component pinned with `constant_other()` cannot drift, so listing it
    # would bury the ones that can. This is the same set that
    # `warn_other_components_fixed()` names.
    expect_length(steady_drift_report(residual_params)$other, 0)
    p <- setComponent(residual_params, "fixed", initial_value = c(1, 2),
                      dynamics_fun = "constant_other")
    expect_length(steady_drift_report(p)$other, 0)

    p <- setComponent(p, "still", initial_value = c(1, 2),
                      dynamics_fun = "settled_component")
    report <- steady_drift_report(p)
    expect_named(report$other, "still")
    expect_equal(report$other[["still"]], 0)
})

test_that("steady_drift_report() measures a component's rate of change", {
    report <- steady_drift_report(component_params)
    expect_named(report$other, "drifter")
    # `exp(2 dt)` per step is a relative rate of 2 per year.
    expect_equal(report$other[["drifter"]], 2, tolerance = 1e-3)
})

test_that("an unmeasurable component state is reported as such", {
    p <- setComponent(residual_params, "opaque",
                      initial_value = list(a = "not numeric"),
                      dynamics_fun = "opaque_component")
    report <- steady_drift_report(p)
    expect_named(report$other, "opaque")
    expect_true(is.na(report$other[["opaque"]]))
    # It is not evidence of motion, so it is not reported as motion. The
    # component's dynamics here return the state unchanged, and a non-numeric
    # state yields `NA` either way, so calling that a drift would call this
    # model unsteady for ever.
    groups <- component_groups(report)
    expect_length(groups$moving, 0)
    expect_identical(groups$unmeasured, "opaque")
    expect_identical(component_drift_txt(report), "")
    expect_match(component_unmeasured_txt(report), "could not be determined")
    expect_match(component_unmeasured_txt(report), "`opaque`")
})

test_that("an unsupported component makes whole-model steadiness unknown", {
    p <- setComponent(residual_params, "opaque",
                      initial_value = list(a = "not numeric"),
                      dynamics_fun = "opaque_component")
    expect_identical(isSteady(p), NA)
    # The report says what could not be checked, and claims nothing beyond it.
    w <- expect_warning(warn_if_not_steady(p, "Context."),
                        "could not be determined")
    expect_no_match(conditionMessage(w), "is not at its steady state")
    expect_no_match(conditionMessage(w), "but the model as a whole is not")
})

test_that("a moving component and an unmeasurable one are reported apart", {
    p <- setComponent(component_params, "opaque",
                      initial_value = list(a = "not numeric"),
                      dynamics_fun = "opaque_component")
    report <- steady_drift_report(p)
    groups <- component_groups(report)
    expect_named(groups$moving, "drifter")
    expect_identical(groups$unmeasured, "opaque")
    w <- expect_warning(warn_if_not_steady(p, "Context."),
                        "not at its steady state")
    expect_match(conditionMessage(w), "component:drifter")
    expect_match(conditionMessage(w), "component:opaque")
})

test_that("the component sentence suits a caller that advances them", {
    report <- steady_drift_report(component_params)
    expect_match(component_drift_txt(report),
                 "steady-state machinery does not settle them")
    # `project_until_settled()` does advance them, so the standing sentence
    # would contradict what the run just did.
    advanced <- component_drift_txt(report, advanced = TRUE)
    expect_match(advanced, "this run does advance them")
    expect_no_match(advanced, "does not settle them")
})

# The exclusion of other components ----

test_that("the biomass drift excludes the other components", {
    # The heart of the change made for #589: a component with dynamics of its
    # own no longer enters the scalar the tolerances are stated against, so a
    # model whose consumers and resource have settled reads as steady even
    # while the component moves at 2 per year.
    expect_equal(steady_biomass_drift(component_params),
                 steady_biomass_drift(residual_params))
    expect_false(isSteady(component_params))
    expect_lt(steady_biomass_drift(component_params), steady_residual_tol())
})

test_that("the excluded component is still measured and reported", {
    # Excluding it from the scalar must not mean not measuring it: the rate is
    # still there for the user to find, and `steady_rates()` still steps the
    # component to get it.
    other <- attr(getSteadyResidual(component_params), "other")
    expect_named(other, "drifter")
    expect_equal(other$drifter, rep(2, 3), tolerance = 1e-3)
    expect_equal(max(abs(other$drifter)),
                 steady_drift_report(component_params)$other[["drifter"]])
})

test_that("a drifting component is assessed when the consumers are steady", {
    warning <- expect_warning(warn_if_not_steady(component_params, "Context."),
                              "not at its steady state")
    expect_match(conditionMessage(warning), "component:drifter/relative_change")
    expect_match(conditionMessage(warning), "getSteadyAssessment")
    expect_no_match(conditionMessage(warning), "Missing convergence support")
})

test_that("a drifting component is named alongside a drifting species", {
    p <- component_params
    initialN(p)[1, ] <- initialN(p)[1, ] * 3
    w <- expect_warning(warn_if_not_steady(p, "Context."),
                        "not at its steady state")
    expect_warning(warn_if_not_steady(p, "Context."), "drifter")
})

test_that("a component with a passing assessment is not named", {
    p <- setComponent(residual_params, "still", initial_value = c(1, 2),
                      dynamics_fun = "settled_component",
                      convergence_fun = "settled_component_assessment")
    expect_identical(component_drift_txt(steady_drift_report(p)), "")
    expect_silent(warn_if_not_steady(p, "Context."))
})

test_that("steady_total_drift() folds the components back in", {
    # The one caller that needs them: `project_until_settled()`, where the
    # components are live rather than pinned.
    report <- steady_drift_report(component_params)
    expect_equal(steady_total_drift(report), 2, tolerance = 1e-3)
    expect_gt(steady_total_drift(report), report$drift)
    # With nothing to fold in it is just the drift.
    expect_identical(steady_total_drift(steady_drift_report(residual_params)),
                     steady_biomass_drift(residual_params))
})

# warn_if_not_steady ----

test_that("warn_if_not_steady() reports unsettled models", {
    expect_warning(warn_if_not_steady(off_steady_params, "Context."),
                   "not at its steady state")
    expect_silent(warn_if_not_steady(residual_params, "Context."))
})

test_that("warn_if_not_steady() honours context and info_level", {
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
    # nolint start: object_name_linter.
    isSteady.DummyModel <- function(params, ...) TRUE
    # nolint end: object_name_linter.
    expect_true(isSteady(dummy))
})

test_that("steady diagnostics share one complete effective rate evaluation", {
    calls <- 0L
    rlang::local_bindings(
        diagnostic_full_rates = function(params, ...) {
            calls <<- calls + 1L
            rates <- mizerRates(params, ...)
            rates$rdd <- rates$rdd / 2
            rates
        },
        diagnostic_fishing_pool = function(n_other, component, rates, dt,
                                           predicted_state, endpoint_rates,
                                           ...) {
            if (is.null(rates$f_mort)) stop("missing fishing mortality")
            stopifnot(is.null(predicted_state), is.null(endpoint_rates))
            n_other[[component]] + dt * sum(rates$f_mort)
        },
        diagnostic_fishing_resource = function(n_pp, rates, dt, predicted_state,
                                               endpoint_rates, ...) {
            if (is.null(rates$f_mort)) stop("missing fishing mortality")
            stopifnot(is.null(predicted_state), is.null(endpoint_rates))
            n_pp * (1 + dt * sum(rates$f_mort))
        }, .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = 2,
                      dynamics_fun = "diagnostic_fishing_pool")
    p@resource_dynamics <- "diagnostic_fishing_resource"
    p@rates_funcs$Rates <- "diagnostic_full_rates"
    dt <- 0.25
    measured <- steady_rates(p, dt = dt)
    expect_identical(calls, 1L)

    rates <- mizerRates(p, n = p@initial_n, n_pp = p@initial_n_pp,
                        n_other = p@initial_n_other, effort = p@initial_effort,
                        rates_fns = projectRateFunctions(p))
    fishing <- sum(rates$f_mort)
    expect_gt(fishing, 0)
    expect_equal(measured$other$pool, fishing / 2)
    present <- p@initial_n_pp > 0
    expect_equal(unname(measured$dn_pp_dt[present] / p@initial_n_pp[present]),
                 rep(fishing, sum(present)), tolerance = 1e-10)
    expect_true(all(measured$dn_pp_dt[!present] == 0))

    # The full-rate provider changes final recruitment. The consumer equation
    # must use that value instead of recomputing the unmodified RDD separately.
    expected <- -consumer_residual(p, n = p@initial_n,
                                   n_pp = p@initial_n_pp,
                                   n_other = p@initial_n_other,
                                   effort = p@initial_effort,
                                   rates = rates, rdd = rates$rdd / 2)
    populated <- p@initial_n > 0
    expect_equal(measured$dNdt[populated], expected[populated])
})

test_that("steady diagnostics use equations without an artificial step", {
    rlang::local_bindings(
        diagnostic_component_rate = function(n_other, component, ...) {
            2 * n_other[[component]]
        },
        diagnostic_unused_step = function(...) stop("step should not be used"),
        .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = c(stock = 2,
                                                                 empty = 0),
                      dynamics_fun = "diagnostic_unused_step",
                      rate_fun = "diagnostic_component_rate")
    expected <- c(stock = 2, empty = NA_real_)
    # The equation is sufficient even if the unused custom step is unavailable
    # in a session that only needs this diagnostic.
    rm("diagnostic_unused_step", envir = .GlobalEnv)
    expect_equal(steady_rates(p, dt = 1)$other$pool, expected)
    expect_equal(steady_rates(p, dt = 1e-8)$other$pool, expected)

    # Tuning temporarily replaces the active step. Its frozen state must not
    # be diagnosed using the original unfrozen derivative from the snapshot.
    p@other_dynamics$pool <- "constant_other"
    expect_equal(steady_rates(p)$other$pool, c(stock = 0, empty = NA_real_))
})

test_that("steady diagnostics reject bad derivatives and propagate errors", {
    rlang::local_bindings(
        diagnostic_bad_rate = function(...) stop("component derivative failed"),
        diagnostic_wrong_shape = function(...) 1,
        .env = .GlobalEnv
    )
    p <- setComponent(NS_params_small, "pool", initial_value = c(a = 1, b = 2),
                      rate_fun = "diagnostic_bad_rate")
    expect_error(steady_rates(p), "component derivative failed")
    p <- setComponent(p, "pool", rate_fun = "diagnostic_wrong_shape")
    expect_error(steady_rates(p), "pool")
})

test_that("Conditional isSteady verdicts retain scope and fixed owners", {
    p <- setComponent(residual_params, "opaque", initial_value = list(a = 1),
                      dynamics_fun = "constant_other")
    verdict <- isSteady(p, scope = "conditional", fixed = "component:opaque")
    expect_true(verdict)
    expect_identical(attr(verdict, "scope"), "conditional")
    expect_identical(attr(verdict, "fixed"), "component:opaque")
})

test_that("Steadiness checks preserve callback errors and assessment time", {
    rlang::local_bindings(
        timed_assessment = function(t, ...) {
            data.frame(measure = "time", value = t,
                       tolerance = 1, units = "year")
        },
        failed_assessment = function(...) stop("assessment programming error"),
        .env = .GlobalEnv
    )
    p <- setComponent(residual_params, "clock", initial_value = 0,
                      dynamics_fun = "constant_other",
                      convergence_fun = "timed_assessment")
    expect_true(isSteady(p, t = 0))
    expect_false(isSteady(p, t = 2))
    expect_warning(warn_if_not_steady(p, "Context.", t = 2), "component:clock")
    p <- setComponent(p, "clock", convergence_fun = "failed_assessment")
    expect_error(isSteady(p), "assessment programming error")
    expect_error(warn_if_not_steady(p, "Context."),
                 "assessment programming error")
    withr::local_options(mizer_info_level = 0)
    expect_error(warn_if_not_steady(p, "Context."),
                 "assessment programming error")
})


test_that("Steady diagnostics identify malformed effective rate fields", {
    rlang::local_bindings(
        nonfinite_assessment_rates = function(params, ...) {
            rates <- mizerRates(params, ...)
            rates$encounter[1, 1] <- NaN
            rates
        },
        nonnumeric_assessment_rates = function(params, ...) {
            rates <- mizerRates(params, ...)
            rates$rdd <- "invalid recruitment"
            rates
        }, .env = .GlobalEnv
    )
    p <- NS_params_small
    p$rates_funcs$Rates <- "nonfinite_assessment_rates"
    expect_error(steady_rates(p),
                 "nonfinite_assessment_rates.*returned non-finite.*encounter")
    expect_error(isSteady(p), "returned non-finite values")
    p$rates_funcs$Rates <- "nonnumeric_assessment_rates"
    expect_error(steady_rates(p),
                 "nonnumeric_assessment_rates.*non-numeric.*rdd")
})
