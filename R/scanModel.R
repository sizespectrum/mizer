# scanModel: scan a model over a range of values for one of its aspects
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' Scan a model over a range of values
#'
#' `r lifecycle::badge("experimental")`
#' Varies one aspect of a model over a range of values and measures a quantity
#' at each of them. At every value the model is projected until it settles onto
#' an attractor, and the quantity is measured on that attractor rather than at
#' whatever state the projection happened to stop at.
#'
#' You say what to vary by giving a function that changes the model, and what to
#' measure by giving a function that computes a quantity from a simulation. So a
#' yield-versus-fishing-mortality curve, a bifurcation diagram over fishing
#' effort and a scan over the resource carrying capacity are all the same
#' function call with different arguments.
#'
#' This is a generic function with a method for objects of class
#' [MizerParams][mizer::MizerParams].
#'
#' @section What is measured, and where:
#'
#' At each scan value the model is projected with [projectUntilSettled()], which
#' stops as soon as it recognises that the model has settled onto a fixed point
#' or onto a limit cycle, and reports which of the two happened. What happens
#' next depends on that answer:
#'
#' \describe{
#'   \item{A fixed point}{The state does not change, so there is nothing to
#'     average. The quantity is read off the settled state with no further
#'     projection at all, and the reported minimum and maximum are equal to it.}
#'   \item{A limit cycle}{The model is projected for **exactly one period** of
#'     the detected cycle and the quantity is averaged over it, which is its
#'     long-term average. The minimum and maximum over the cycle are reported
#'     too.}
#'   \item{Neither}{The model did not settle within `t_max` years, or a species
#'     went extinct. The quantity is averaged over the last `t_sample` years and
#'     the scan values concerned are named in a message, because those points
#'     should not be relied on.}
#' }
#'
#' Averaging over exactly one period is both faster and more accurate than
#' averaging over a fixed number of years. A window that is not a whole number
#' of periods leaves a residue of the oscillation in the average, which shows up
#' as a jagged curve. The window is rounded to a whole number of time steps, so
#' it can differ from the true period by up to `dt/2`; if you need the average
#' more accurately, reduce `dt` rather than lengthening the window, because a
#' longer window that is not a whole number of periods is worse, not better.
#'
#' @section Writing the two functions:
#'
#' `set_func(params, value)` takes a `MizerParams` object and one entry of
#' `scan_values` and returns a modified `MizerParams`. It must be **idempotent**
#' — `set_func(set_func(p, v), v)` must give the same thing as
#' `set_func(p, v)` — because with `continuation = TRUE` it is applied to the
#' object it returned at the previous scan value. Setting something is
#' idempotent; appending something is not, so a function that adds a gear must
#' check whether the gear is already there. See [scanFishingMortality()] for a
#' worked example.
#'
#' There is no `effort` argument, because there does not need to be one:
#' [project()] and [projectUntilSettled()] both take the fishing effort from
#' `params@initial_effort`, so a `set_func()` that changes the effort is all it
#' takes to scan over effort, and a scan over something else never has to
#' mention fishing at all.
#'
#' `value_func(sim)` takes a `MizerSim` and returns either a time by series
#' matrix, as [getBiomass()], [getYield()], [getSSB()], [getN()] and
#' [sizeIntegral()] all do, or a plain numeric vector over time, as
#' [getMeanWeight()] does. When it returns a matrix carrying `value_name` and
#' `units` attributes — which all of mizer's `MizerSim` methods do — those are
#' used for the y-axis label unless you override them.
#'
#' Note that at a fixed point `value_func()` is handed a simulation with a
#' single time step, so a function that needs more than one time step will not
#' work there. Set `sample_all = TRUE` to force the sampling projection at every
#' scan value.
#'
#' Neither function can be given extra arguments through `...`, which is
#' reserved for `distance_func`. Use a closure instead, for example
#' `value_func = function(sim) getBiomass(sim, min_w = 10)`.
#'
#' @param params An object of class `MizerParams`.
#' @param scan_values A numeric vector of values to scan over.
#' @param set_func A function of `(params, value)` returning a modified
#'   `MizerParams` object. Several are provided: [scanEffort()],
#'   [scanFishingMortality()] and [scanSpeciesParam()].
#' @param value_func A function of a `MizerSim` returning the quantity to
#'   measure. Defaults to [getBiomass()].
#' @param species The species to keep in the result. By default all of the
#'   series that `value_func` returns.
#' @param scan_name A string naming the quantity that is varied, used for the
#'   x-axis label and as the name of the first column of the result. Taken from
#'   `set_func` when it supplies one.
#' @param scan_units A string giving the units of that quantity.
#' @param value_name A string naming the quantity that is measured, used for the
#'   y-axis label and as the name of the second column of the result. Taken from
#'   what `value_func` returns when it supplies one.
#' @param value_units A string giving the units of that quantity.
#' @param reference_lines An optional named numeric vector of positions on the x
#'   axis for [plot.MizerScan()] to mark with vertical dashed lines, for example
#'   `c(F_MSY = 0.32)`. Taken from `set_func` when it supplies one.
#' @param current_scan_value The value at which the model currently sits. When
#'   given, the scan works outwards from it in both directions so that every
#'   projection starts from a neighbouring attractor rather than from a distant
#'   state, and each of the two directions begins again at the model as it was
#'   given. Pass `"auto"` to ask `set_func` for it. By default the values are
#'   scanned in the order given, which is also how you trace a hysteresis branch
#'   deliberately: pass a decreasing `scan_values`.
#' @param continuation Whether each scan value should start from the attractor
#'   reached at the previous one. Default TRUE.
#' @param distance_func A function that will be called after every `t_per` years
#'   with both the previous and the new state and that should return a number
#'   measuring the distance between them. See [distanceSSLogN()].
#' @param tol The projection at each scan value stops once the number returned
#'   by `distance_func` for two states `t_per` years apart drops below `tol`.
#'   The default is tighter than the one [projectUntilSettled()] uses on its own,
#'   because a scan produces a curve, and a loosely converged point does not
#'   average away: it shows up as a kink in the curve and as spurious width in
#'   the band. Loosen it to go faster, and use the `residual` column of the
#'   result to check whether that cost you anything.
#' @param t_per The interval in years at which convergence is checked. Should be
#'   an odd multiple of `dt`.
#' @param t_max The longest time to project at each scan value.
#' @param dt The time step to use.
#' @param t_save The interval at which the biomass summary used for limit-cycle
#'   detection is recorded, see [projectUntilSettled()].
#' @param amplitude_tol The minimum relative biomass amplitude for a persistent
#'   oscillation to count as a limit cycle rather than a fixed point.
#' @param amp_rel_tol Maximum relative change of amplitude between successive
#'   periods for a cycle to count as settled.
#' @param extinction_threshold A species is treated as going extinct once its
#'   reproduction rate falls below this fraction of its value at the start of
#'   the projection.
#' @param method The numerical method to use, see [project()].
#' @param t_sample The number of years over which to average when the model has
#'   settled onto neither a fixed point nor a limit cycle.
#' @param sample_all Whether to run the sampling projection even at a fixed
#'   point, where it is not otherwise needed. Set this if `value_func` needs
#'   more than one time step.
#' @param progress_bar Whether to show a text progress bar over the scan values.
#' @param info_level Controls how much the projections say for themselves.
#'   Defaults to 0, because a scan makes one projection per scan value and
#'   summarises them itself; raise it when investigating why one of them
#'   behaved oddly.
#' @param ... Further arguments are passed on to `distance_func`.
#'
#' @return An object of class [MizerScan], which is a data frame with one row
#'   per scan value and series, carrying the metadata that [plot()] needs.
#'
#' @family scan functions
#' @concept scan
#' @seealso [MizerScan()], [plot.MizerScan()], [scanEffort()],
#'   [scanFishingMortality()], [scanSpeciesParam()], [plotYieldVsF()]
#' @export
#' @examples
#' \donttest{
#' # A bifurcation diagram over fishing effort
#' scan <- scanModel(NS_params, scan_values = seq(0, 2, 0.25),
#'                   set_func = scanEffort(), value_func = getYield)
#' plot(scan, style = "envelope")
#'
#' # A yield curve for a single species, and the F at which it is largest
#' cod <- scanModel(NS_params, scan_values = seq(0, 1.2, 0.2),
#'                  set_func = scanFishingMortality("Cod"),
#'                  value_func = getYield, species = "Cod")
#' plot(cod, mark_max = TRUE, log_y = FALSE)
#' attr(cod, "at_max")
#'
#' # Scanning something that has nothing to do with fishing
#' kappa <- resource_params(NS_params)$kappa
#' plot(scanModel(NS_params, scan_values = kappa * c(0.5, 1, 2),
#'                set_func = function(params, value) {
#'                    resource_params(params)$kappa <- value
#'                    params
#'                },
#'                scan_name = "Resource capacity", scan_units = "g"),
#'      log_x = TRUE)
#' }
scanModel <- function(params, scan_values, set_func,
                      value_func = getBiomass,
                      species = NULL,
                      scan_name = NULL, scan_units = NULL,
                      value_name = NULL, value_units = NULL,
                      reference_lines = NULL,
                      current_scan_value = NULL,
                      continuation = TRUE,
                      distance_func = distanceSSLogN,
                      tol = 0.001, t_per = 1.5, t_max = 100, dt = 0.1,
                      t_save = dt,
                      amplitude_tol = 0.01, amp_rel_tol = 0.1,
                      extinction_threshold = 1e-6,
                      method = c("euler", "predictor_corrector", "tr_bdf2"),
                      t_sample = 10, sample_all = FALSE,
                      progress_bar = interactive(),
                      info_level = 0, ...) {
    UseMethod("scanModel")
}

#' @export
scanModel.MizerParams <- function(params, scan_values, set_func,
                                  value_func = getBiomass,
                                  species = NULL,
                                  scan_name = NULL, scan_units = NULL,
                                  value_name = NULL, value_units = NULL,
                                  reference_lines = NULL,
                                  current_scan_value = NULL,
                                  continuation = TRUE,
                                  distance_func = distanceSSLogN,
                                  tol = 0.001, t_per = 1.5, t_max = 100,
                                  dt = 0.1, t_save = dt,
                                  amplitude_tol = 0.01, amp_rel_tol = 0.1,
                                  extinction_threshold = 1e-6,
                                  method = c("euler", "predictor_corrector",
                                             "tr_bdf2"),
                                  t_sample = 10, sample_all = FALSE,
                                  progress_bar = interactive(),
                                  info_level = 0, ...) {
    params <- validParams(params)
    method <- match.arg(method)
    assert_that(is.numeric(scan_values),
                is.function(set_func),
                is.function(value_func),
                is.number(tol), tol > 0,
                is.number(t_per), t_per > 0,
                is.number(t_max), t_max > 0,
                is.number(dt), dt > 0,
                is.number(t_sample), t_sample > 0,
                is.number(amplitude_tol), amplitude_tol > 0,
                is.number(extinction_threshold), extinction_threshold >= 0,
                is.flag(continuation),
                is.flag(sample_all),
                is.flag(progress_bar))
    if (length(scan_values) == 0) {
        stop("`scan_values` must contain at least one value.")
    }
    if (anyNA(scan_values)) {
        stop("`scan_values` must not contain NA.")
    }

    # Labels and reference lines can be supplied by the setter, so that
    # scanEffort() and friends produce a properly labelled plot on their own.
    scan_name <- scan_name %||% attr(set_func, "scan_name") %||% "Scan value"
    scan_units <- scan_units %||% attr(set_func, "scan_units")
    if (is.null(reference_lines)) {
        ref_func <- attr(set_func, "reference_lines")
        if (is.function(ref_func)) reference_lines <- ref_func(params)
    }

    if (identical(current_scan_value, "auto")) {
        current_func <- attr(set_func, "current_scan_value")
        if (!is.function(current_func)) {
            stop('current_scan_value = "auto" needs a `set_func` that says ',
                 "where the model currently sits, which this one does not. ",
                 "Give the value explicitly instead.")
        }
        current_scan_value <- current_func(params)
    }
    if (!is.null(current_scan_value)) {
        assert_that(is.number(current_scan_value))
        if (!continuation) {
            signal_info("current_scan_value",
                        paste0("`current_scan_value` only affects the order ",
                               "in which the scan values are projected, which ",
                               "does not matter when `continuation = FALSE`. ",
                               "It will be ignored."),
                        severity = "warning", unhandled = "show")
        }
    }

    arms <- sweep_arms(scan_values, current_scan_value)

    pb <- NULL
    if (isTRUE(progress_bar)) {
        pb <- utils::txtProgressBar(min = 0, max = length(scan_values),
                                    style = 3)
        on.exit(close(pb), add = TRUE)
    }

    rows <- vector("list", length(scan_values))
    types <- rep(NA_character_, length(scan_values))
    series_names <- NULL
    step <- 0L

    for (arm in arms) {
        # Each arm starts again from the model as it was given, so that the
        # arm working upwards from the current value does not begin at the
        # attractor reached at the far end of the arm working downwards.
        p_run <- params
        for (i in arm) {
            step <- step + 1L
            p <- set_func(p_run, scan_values[[i]])
            if (!is(p, "MizerParams")) {
                stop("`set_func` must return a MizerParams object, but at the ",
                     "scan value ", signif(scan_values[[i]], 3),
                     " it returned an object of class ",
                     paste(class(p), collapse = "/"), ".")
            }

            settled <- project_until_settled(
                p, distance_func = distance_func,
                t_per = t_per, t_max = t_max, dt = dt, t_save = t_save,
                tol = tol, amplitude_tol = amplitude_tol,
                amp_rel_tol = amp_rel_tol,
                extinction_threshold = extinction_threshold,
                method = method, return_sim = FALSE, progress_bar = FALSE,
                info_level = info_level, ...
            )
            conv <- attr(settled, "convergence")

            measured <- measure_on_attractor(
                settled, value_func, conv, dt = dt, t_sample = t_sample,
                sample_all = sample_all, method = method,
                default_name = value_name %||% "Value"
            )
            if (is.null(series_names)) {
                series_names <- measured$series
                value_name <- value_name %||% measured$value_name
                value_units <- value_units %||% measured$value_units
                value_type <- measured$type
                # Checked as soon as the series are known, so that a
                # mistyped name is not paid for with the whole scan first.
                if (!is.null(species)) select_scan_series(series_names, species)
            } else if (!identical(measured$series, series_names)) {
                stop("`value_func` returned different series at different ",
                     "scan values: ", paste(series_names, collapse = ", "),
                     " and then ", paste(measured$series, collapse = ", "),
                     ". A scan needs the same series throughout.")
            }

            rows[[i]] <- data.frame(
                x = scan_values[[i]],
                value = measured$mean,
                Species = series_names,
                ymin = measured$min,
                ymax = measured$max,
                type = conv$type,
                settled = isTRUE(conv$settled),
                period = conv$period %||% NA_real_,
                residual = conv$residual %||% NA_real_,
                row.names = NULL,
                stringsAsFactors = FALSE
            )
            types[[i]] <- conv$type %||% NA_character_

            if (continuation) p_run <- settled
            if (!is.null(pb)) utils::setTxtProgressBar(pb, step)
        }
    }

    frame <- do.call(rbind, rows)
    names(frame)[1:2] <- c(scan_name, value_name %||% "Value")

    if (!is.null(species)) {
        # Selected against the series the value function actually returned,
        # which need not be species at all.
        frame <- frame[select_scan_series(frame$Species, species), ,
                       drop = FALSE]
        rownames(frame) <- NULL
    }

    report_scan_convergence(scan_values, types, scan_name, t_max, t_sample)

    MizerScan(frame,
              scan_name = scan_name, scan_units = scan_units,
              value_name = value_name %||% "Value", value_units = value_units,
              type = value_type, params = params,
              reference_lines = reference_lines,
              settings = list(tol = tol, t_per = t_per, t_max = t_max,
                              dt = dt, t_sample = t_sample,
                              continuation = continuation,
                              method = method))
}

#' Measure a quantity on the attractor a projection settled on
#'
#' @param settled The MizerParams returned by `project_until_settled()`.
#' @param value_func The function measuring the quantity.
#' @param conv The `"convergence"` attribute of `settled`.
#' @param dt The time step.
#' @param t_sample The averaging window to use when nothing settled.
#' @param sample_all Whether to sample even at a fixed point.
#' @param method The numerical method.
#' @param default_name The series name to use when `value_func` supplies none.
#' @return A list with the mean, minimum and maximum over the attractor, the
#'   names of the series and the metadata read off `value_func`'s result.
#' @keywords internal
measure_on_attractor <- function(settled, value_func, conv, dt, t_sample,
                                 sample_all, method, default_name = "Value") {
    whole_period <- identical(conv$type, "cycle") &&
        is.finite(conv$period) && conv$period > 0
    fixed_point <- identical(conv$type, "below_tolerance") && !sample_all

    if (fixed_point) {
        # Nothing moves, so a snapshot carries all the information a longer
        # run would, and costs no projection at all.
        vals <- value_func(params_as_sim(settled))
        keep <- 1L
    } else {
        window <- if (whole_period) conv$period else t_sample
        # project() advances in whole time steps.
        no_steps <- max(1, round(window / dt))
        sim <- project(settled, t_max = no_steps * dt, dt = dt, t_save = dt,
                       method = method, progress_bar = FALSE)
        vals <- value_func(sim)
        # With t_save = dt there are no_steps + 1 saved times, 0 to
        # no_steps * dt. Over exactly one period the last of them repeats the
        # phase of the first, so it must be left out of the mean. It stays in
        # the minimum and the maximum, which are over the whole window.
        n_rows <- if (is.null(dim(vals))) length(vals) else nrow(vals)
        keep <- if (whole_period && n_rows > 1) {
            seq_len(n_rows - 1L)
        } else {
            seq_len(n_rows)
        }
    }

    mat <- as_series_matrix(vals, default_name = default_name)
    list(series = colnames(mat),
         mean = colMeans(mat[keep, , drop = FALSE]),
         min = apply(mat, 2, min),
         max = apply(mat, 2, max),
         value_name = attr(vals, "value_name"),
         value_units = attr(vals, "units"),
         type = resolve_array_type(attr(vals, "type"),
                                   attr(vals, "value_name"),
                                   attr(vals, "units")))
}

#' Bring what a value function returned into a time by series matrix
#'
#' @param x The return value of a `value_func`.
#' @param default_name The column name to use for a single unnamed series.
#' @return A matrix with time in the rows and the series in the columns.
#' @keywords internal
as_series_matrix <- function(x, default_name = "Value") {
    if (is.data.frame(x)) {
        stop("`value_func` returned a data frame. A scan needs one value per ",
             "series per time, so pick out the column you want, for example\n",
             '  value_func = function(sim) getCommunitySlope(sim)[, "slope"]')
    }
    if (is.null(dim(x))) {
        if (!is.numeric(x)) {
            stop("`value_func` must return numbers, but it returned an object ",
                 "of class ", paste(class(x), collapse = "/"), ".")
        }
        return(matrix(x, ncol = 1,
                      dimnames = list(names(x), default_name)))
    }
    if (length(dim(x)) != 2) {
        stop("`value_func` must return a matrix with time in the rows and the ",
             "series in the columns, or a numeric vector over time, but it ",
             "returned an array with ", length(dim(x)), " dimensions.")
    }
    if (is.null(colnames(x))) {
        colnames(x) <- if (ncol(x) == 1) {
            default_name
        } else {
            paste("Series", seq_len(ncol(x)))
        }
    }
    x
}

#' The order in which to project the scan values
#'
#' With continuation each scan value starts from the attractor reached at the
#' previous one, so it pays to visit them in an order where consecutive values
#' are close together. When the value the model currently sits at is known, that
#' means working outwards from it in both directions; otherwise the order the
#' user gave is used, which is also what lets a decreasing `scan_values` trace a
#' hysteresis branch deliberately.
#'
#' The two directions are returned as separate arms rather than as one
#' sequence, because each has to begin again at the model as it was given. Run
#' as one sequence they would carry the attractor from the far end of the
#' descending arm into the start of the ascending arm, which is the opposite of
#' starting each projection from a neighbour, and in a model with coexisting
#' attractors would follow the wrong branch.
#'
#' @param scan_values The values to scan over.
#' @param current_scan_value The value the model currently sits at, or NULL.
#' @return A list of integer vectors, each holding indices into `scan_values` in
#'   the order they should be projected. Each arm is to be started from the
#'   unmodified model.
#' @keywords internal
sweep_arms <- function(scan_values, current_scan_value = NULL) {
    if (is.null(current_scan_value)) {
        return(list(seq_along(scan_values)))
    }
    below <- which(scan_values < current_scan_value)
    above <- which(scan_values >= current_scan_value)
    arms <- list(rev(below[order(scan_values[below])]),
                 above[order(scan_values[above])])
    arms[lengths(arms) > 0]
}

#' Report the scan values that did not settle on a fixed point
#'
#' @param scan_values The values that were scanned.
#' @param types The attractor type reached at each value, one per value.
#' @param scan_name The name of the scanned quantity.
#' @param t_max The time limit that was used.
#' @param t_sample The averaging window that was used.
#' @return Nothing; called for its messages.
#' @keywords internal
report_scan_convergence <- function(scan_values, types, scan_name, t_max,
                                    t_sample) {
    # `which()` rather than the logical vector itself, so that a missing type
    # cannot put an NA among the scan values that are named.
    name_them <- function(type) {
        paste(signif(scan_values[which(types == type)], 3), collapse = ", ")
    }
    # `...` is only forced inside the `if`, so a report that is not needed
    # costs nothing to have written down here.
    report <- function(type, ...) {
        if (any(types == type, na.rm = TRUE)) {
            signal_info("scan", paste0(...), unhandled = "show")
        }
    }
    report("cycle",
           "The model settled onto a limit cycle at ", scan_name, " = ",
           name_them("cycle"),
           ". The value there is the average over one period of the cycle.")
    report("not_converged",
           "The model did not settle onto an attractor within ", t_max,
           " years at ", scan_name, " = ", name_them("not_converged"),
           ". The value there is only the average over the last ", t_sample,
           " years and should not be relied on. Increase `t_max` to give the ",
           "model more time to settle, or loosen `tol`. The `residual` ",
           "column says how fast the abundances were still changing, in ",
           "1/year.")
    report("extinction",
           "A species was on its way to extinction at ", scan_name, " = ",
           name_them("extinction"),
           ", which stopped the projection early. The value there is only ",
           "the average over the ", t_sample, " years following that point.")
    invisible(NULL)
}
