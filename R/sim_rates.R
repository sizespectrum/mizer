# Machinery for evaluating rates over the saved time steps of a MizerSim.
#
# These helpers are not rate functions themselves: they select saved time
# steps, assemble the model state at a saved step, ask for just the rates a
# caller needs, and wrap the result in the appropriate array class. They are
# shared infrastructure, used by the `MizerSim` rate methods in
# rate_functions.R but also by summary_methods.R, diffusion.R, steadyNewton.R,
# getRequiredRDD.R, plots.R, animateSpectra.R, MizerSim-class.R and
# setInitialValues.R.

#' Get selected saved time steps for a simulation rate
#'
#' Internal helper used by `MizerSim` rate methods. If `time_range` is missing,
#' all saved simulation times are selected; otherwise the request is delegated to
#' [get_time_elements()].
#'
#' @param sim A `MizerSim` object.
#' @param time_range A numeric or character vector of times.
#'
#' @return A named logical vector indicating the selected saved time steps.
#' @keywords internal
get_sim_rate_time_elements <- function(sim, time_range) {
    # A missing time range means "all saved times", matching the public
    # MizerSim rate methods.
    if (missing(time_range)) {
        time_range <- dimnames(sim@n)$time
    }
    get_time_elements(sim, time_range)
}

#' Extract one saved simulation state for a rate calculation
#'
#' Internal helper used by `MizerSim` rate methods to rebuild the single-time
#' inputs expected by `MizerParams` rate methods.
#'
#' @param sim A `MizerSim` object.
#' @param time_idx Integer index of the saved time step to extract.
#'
#' @return A list with entries `n`, `n_pp`, `n_other`, `effort`, and `t`.
#' @keywords internal
get_sim_rate_slice <- function(sim, time_idx) {
    # Rebuild `n` explicitly so that a one-species simulation still gives a
    # species x size matrix rather than a vector.
    n <- array(sim@n[time_idx, , ], dim = dim(sim@n)[2:3])
    dimnames(n) <- dimnames(sim@n)[2:3]

    # The `n_other` slot is a time x component list-array. A single row needs
    # its component names restored before passing it to rate functions.
    n_other <- sim@n_other[time_idx, ]
    names(n_other) <- dimnames(sim@n_other)$component

    list(
        n = n,
        n_pp = sim@n_pp[time_idx, ],
        n_other = n_other,
        effort = sim@effort[time_idx, ],
        t = as.numeric(dimnames(sim@n)$time[[time_idx]])
    )
}

#' Apply a species-by-size rate function over saved simulation times
#'
#' Internal helper used by `MizerSim` rate methods whose one-time result is an
#' `ArraySpeciesBySize`. The helper applies the supplied rate function to each
#' selected time slice, stacks the results, and restores the appropriate mizer
#' array class when dimensions have not been dropped.
#'
#' @param sim A `MizerSim` object.
#' @param time_range A numeric or character vector of times.
#' @param drop If `TRUE`, dimensions of length 1 are dropped from the result.
#' @param rate_fun A function accepting a single simulation slice as returned by
#'   `get_sim_rate_slice()`.
#' @param value_name Name of the value stored in the returned array.
#' @param units Optional units of the value stored in the returned array.
#' @param type The kind of quantity the values are, see [ArraySpeciesBySize()]
#'   and [array_types].
#'
#' @return A time x species x size array, possibly with dimensions dropped.
#' @keywords internal
get_species_size_rate_from_sim <- function(sim, time_range, drop,
                                           rate_fun, value_name,
                                           units = NULL,
                                           type = NULL,
                                           representation = "point") {
    time_elements <- get_sim_rate_time_elements(sim, time_range)

    # Apply the one-time rate calculation to each selected saved time. The
    # result from each call is species x size, so `aaply()` stacks them into a
    # time x species x size array.
    rate_time <- plyr::aaply(which(time_elements), 1, function(time_idx) {
        rate_fun(get_sim_rate_slice(sim, time_idx))
    }, .drop = FALSE)
    names(dimnames(rate_time))[[1]] <- "time"

    result <- rate_time[, , , drop = drop]

    # Restore the richer array classes when the requested dropping has left the
    # dimensions in a shape those classes represent. With a single species and
    # `drop = TRUE`, the result is time x size and should stay a plain matrix.
    if (is.array(result) && length(dim(result)) == 3) {
        result <- ArrayTimeBySpeciesBySize(result,
                                          value_name = value_name,
                                          units = units,
                                          type = type,
                                          params = sim@params,
                                          representation = representation)
    } else if (is.matrix(result) &&
               names(dimnames(result))[[1]] == "sp") {
        result <- ArraySpeciesBySize(result, value_name = value_name,
                                     units = units, type = type,
                                     params = sim@params,
                                     representation = representation)
    }
    result
}

#' Apply a species rate function over saved simulation times
#'
#' Internal helper used by `MizerSim` rate methods whose one-time result is a
#' named vector with one value for each species.
#'
#' @param sim A `MizerSim` object.
#' @param time_range A numeric or character vector of times.
#' @param rate_fun A function accepting a single simulation slice as returned by
#'   `get_sim_rate_slice()`.
#' @param value_name Name of the value stored in the returned array.
#' @param units Optional units of the value stored in the returned array.
#'
#' @return An `ArrayTimeBySpecies` object with dimensions time x species.
#' @keywords internal
get_species_time_rate_from_sim <- function(sim, time_range, rate_fun,
                                           value_name, units = NULL) {
    time_elements <- get_sim_rate_time_elements(sim, time_range)
    time_idx <- which(time_elements)
    species <- sim@params@species_params$species

    # Each one-time rate is a species vector. `vapply()` stacks these as
    # species x time, so transpose to the ArrayTimeBySpecies convention.
    rate_time <- t(vapply(time_idx, function(idx) {
        rate_fun(get_sim_rate_slice(sim, idx))
    }, numeric(length(species))))
    dimnames(rate_time) <- list(time = names(time_elements)[time_elements],
                                sp = species)
    ArrayTimeBySpecies(rate_time, value_name = value_name, units = units,
                       params = sim@params)
}


# Dependency graph of the rate functions: for each rate, the other rates that
# must be calculated before it. The names match the entries of
# `params@rates_funcs`. Used by `mizer_rates_subset()` to work out the minimal
# set of rates a `MizerSim` getter needs at each saved time step.
.rate_dependencies <- list(
    Encounter       = character(0),
    FeedingLevel    = "Encounter",
    EReproAndGrowth = c("Encounter", "FeedingLevel"),
    ERepro          = "EReproAndGrowth",
    EGrowth         = c("ERepro", "EReproAndGrowth"),
    Diffusion       = "FeedingLevel",
    PredRate        = "FeedingLevel",
    PredMort        = "PredRate",
    FMort           = c("EGrowth", "PredMort"),
    Mort            = c("FMort", "PredMort"),
    RDI             = c("ERepro", "EGrowth", "Diffusion", "Mort"),
    RDD             = "RDI",
    ResourceMort    = "PredRate"
)

#' Determine which rates must be calculated to obtain a set of target rates
#'
#' Internal helper returning the transitive closure of `targets` over
#' `.rate_dependencies`, in an order in which the rates can be calculated (each
#' rate appears after all the rates it depends on).
#'
#' @param targets Character vector of rate names (as in `params@rates_funcs`).
#' @return A character vector of rate names.
#' @keywords internal
needed_rates <- function(targets) {
    need <- character(0)
    visit <- function(x) {
        if (x %in% need) return(invisible())
        for (dep in .rate_dependencies[[x]]) visit(dep)
        need[[length(need) + 1L]] <<- x
    }
    for (target in targets) visit(target)
    need
}

#' Calculate a selected subset of the rates
#'
#' Internal helper used by the `MizerSim` rate getters. Given rate functions
#' already resolved once with `projectRateFunctions()`, it calculates only those
#' rates needed to obtain the requested `targets` (plus their dependencies),
#' avoiding both the per-time-step cost of re-resolving the functions and the
#' cost of computing rates that are not required. The individual calculations
#' mirror those in [mizerRates()].
#'
#' @param params A valid `MizerParams` object.
#' @param n A matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size.
#' @param n_other A named list of the abundances of other components.
#' @param t The time for the calculation.
#' @param effort The fishing effort. Only used when a target requires the
#'   fishing mortality.
#' @param rates_fns Named list of resolved rate functions, as returned by
#'   `projectRateFunctions()`.
#' @param targets Character vector of rate names (as in `params@rates_funcs`)
#'   to calculate.
#' @param ... Passed on to the individual rate functions.
#' @return A named list of the calculated rates, using the same element names
#'   as the list returned by [mizerRates()].
#' @keywords internal
mizer_rates_subset <- function(params, n, n_pp, n_other, t, effort,
                               rates_fns, targets, ...) {
    need <- needed_rates(targets)
    r <- list()
    if ("Encounter" %in% need)
        r$encounter <- rates_fns$Encounter(
            params, n = n, n_pp = n_pp, n_other = n_other, t = t, ...)
    if ("FeedingLevel" %in% need)
        r$feeding_level <- rates_fns$FeedingLevel(
            params, n = n, n_pp = n_pp, n_other = n_other,
            encounter = r$encounter, t = t, ...)
    if ("EReproAndGrowth" %in% need)
        r$e <- rates_fns$EReproAndGrowth(
            params, n = n, n_pp = n_pp, n_other = n_other,
            encounter = r$encounter, feeding_level = r$feeding_level,
            t = t, ...)
    if ("ERepro" %in% need)
        r$e_repro <- rates_fns$ERepro(
            params, n = n, n_pp = n_pp, n_other = n_other,
            e = r$e, t = t, ...)
    if ("EGrowth" %in% need)
        r$e_growth <- rates_fns$EGrowth(
            params, n = n, n_pp = n_pp, n_other = n_other,
            e_repro = r$e_repro, e = r$e, t = t, ...)
    if ("Diffusion" %in% need)
        r$diffusion <- rates_fns$Diffusion(
            params, n = n, n_pp = n_pp, n_other = n_other,
            feeding_level = r$feeding_level, t = t, ...)
    if ("PredRate" %in% need)
        r$pred_rate <- rates_fns$PredRate(
            params, n = n, n_pp = n_pp, n_other = n_other,
            feeding_level = r$feeding_level, t = t, ...)
    if ("PredMort" %in% need)
        r$pred_mort <- rates_fns$PredMort(
            params, n = n, n_pp = n_pp, n_other = n_other,
            pred_rate = r$pred_rate, t = t, ...)
    if ("FMort" %in% need)
        r$f_mort <- rates_fns$FMort(
            params, n = n, n_pp = n_pp, n_other = n_other,
            effort = effort, t = t, e_growth = r$e_growth,
            pred_mort = r$pred_mort, ...)
    if ("Mort" %in% need)
        r$mort <- rates_fns$Mort(
            params, n = n, n_pp = n_pp, n_other = n_other,
            f_mort = r$f_mort, pred_mort = r$pred_mort, t = t, ...)
    if ("RDI" %in% need)
        r$rdi <- rates_fns$RDI(
            params, n = n, n_pp = n_pp, n_other = n_other,
            e_growth = r$e_growth, mort = r$mort, diffusion = r$diffusion,
            e_repro = r$e_repro, t = t, ...)
    if ("RDD" %in% need)
        r$rdd <- rates_fns$RDD(
            params, rdi = r$rdi, species_params = params@species_params,
            t = t, ...)
    if ("ResourceMort" %in% need)
        r$resource_mort <- rates_fns$ResourceMort(
            params, n = n, n_pp = n_pp, n_other = n_other,
            pred_rate = r$pred_rate, t = t, ...)
    r
}

#' Build a `MizerSim` rate getter that resolves the rate functions once
#'
#' Internal helper capturing the pattern shared by the `MizerSim` rate getters
#' that return a species-by-size array. It validates the params and resolves the
#' rate functions a single time, then for each saved time step calculates only
#' the required `target` rate with [mizer_rates_subset()] and extracts the
#' element named `slot` from the result.
#'
#' @param sim A `MizerSim` object.
#' @param time_range Passed to the sim iteration helper.
#' @param drop Passed to the sim iteration helper.
#' @param target Name of the rate to calculate (as in `params@rates_funcs`).
#' @param slot Name of the element to extract from the [mizer_rates_subset()]
#'   result (e.g. `"e_growth"` for the `EGrowth` rate).
#' @param value_name,units Metadata for the returned array.
#' @param use_sim_effort If `TRUE`, the saved effort at each time step is used;
#'   otherwise the initial effort is used (matching the behaviour of the
#'   corresponding `MizerParams` getter).
#' @param ... Passed on to the rate functions.
#' @return An `ArrayTimeBySpeciesBySize` object (or a reduced array if `drop`).
#' @keywords internal
sim_size_rate <- function(sim, time_range, drop, target, slot,
                          value_name, units = NULL, type = NULL,
                          use_sim_effort = FALSE,
                          representation = "point", ...) {
    params <- validParams(sim@params)
    rates_fns <- projectRateFunctions(params)
    effort <- params@initial_effort
    get_species_size_rate_from_sim(
        sim, time_range, drop,
        function(slice) {
            r <- mizer_rates_subset(
                params, n = slice$n, n_pp = slice$n_pp,
                n_other = slice$n_other, t = slice$t,
                effort = if (use_sim_effort) slice$effort else effort,
                rates_fns = rates_fns, targets = target, ...)
            m <- r[[slot]]
            # Normalise dimnames to match the corresponding `MizerParams`
            # getter: rates over the community size grid take the `metab`
            # dimnames (as `ArraySpeciesBySize()` would set), while the
            # predation rate runs over the full prey size grid.
            if (identical(dim(m), dim(params@metab))) {
                dimnames(m) <- dimnames(params@metab)
            } else {
                dimnames(m) <- list(
                    sp = dimnames(params@initial_n)$sp,
                    w_prey = as.character(signif(params@w_full, 3)))
            }
            m
        },
        value_name = value_name, units = units, type = type,
        representation = representation)
}

#' Build a `MizerSim` rate getter that resolves the rate functions once
#'
#' Like [sim_size_rate()] but for getters that return one value per species at
#' each time step (a time-by-species array), such as [getRDI()] and [getRDD()].
#' By default these use the initial effort, matching their `MizerParams` counterparts.
#'
#' @inheritParams sim_size_rate
#' @return An `ArrayTimeBySpecies` object.
#' @keywords internal
sim_species_rate <- function(sim, time_range, target, slot,
                             value_name, units = NULL,
                             use_sim_effort = FALSE, ...) {
    params <- validParams(sim@params)
    rates_fns <- projectRateFunctions(params)
    effort <- params@initial_effort
    get_species_time_rate_from_sim(
        sim, time_range,
        function(slice) {
            r <- mizer_rates_subset(
                params, n = slice$n, n_pp = slice$n_pp,
                n_other = slice$n_other, t = slice$t,
                effort = if (use_sim_effort) slice$effort else effort,
                rates_fns = rates_fns, targets = target, ...)
            r[[slot]]
        },
        value_name = value_name, units = units)
}



#' Get array indices for a time range in a MizerSim object
#'
#' Internal helper to select the saved time points whose times lie between the
#' smallest and largest values in `time_range`, inclusive.
#' @param sim A MizerSim object.
#' @param time_range A numeric or character vector of times. Only the range of
#'   values matters, so all saved times between `min(time_range)` and
#'   `max(time_range)` are selected.
#' @param slot_name Obsolete, kept only for backward compatibility with early
#'   versions where different time-based slots could have different time grids.
#'   Leave at the default.
#' @return A named logical vector, with one entry for each saved time in `sim`,
#'   indicating whether that time lies in the requested range.
#' @export
#' @concept helper
get_time_elements <- function(sim, time_range, slot_name = "n") {
    assert_that(is(sim, "MizerSim"))
    time_range <- range(as.numeric(time_range))
    # Check that time range is even in object
    sim_times <- as.numeric(dimnames(slot(sim, slot_name))$time)
    sim_time_range <- range(sim_times)
    if ((time_range[1] < sim_time_range[1]) ||
        (time_range[2] > sim_time_range[2])) {
        stop("Time range is outside the time range of the model.")
    }
    time_elements <- (sim_times >= time_range[1]) & (sim_times <= time_range[2])
    if (!any(time_elements)) {
        stop("The time range does not contain any simulation results.")
    }
    names(time_elements) <- dimnames(sim@effort)$time
    return(time_elements)
}
