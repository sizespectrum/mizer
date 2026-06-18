# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Distributed under the GPL 3 or later
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Note on calling signatures of rate functions
#
# Most rate functions have the arguments as getRates():
#    params, n, n_pp, n_other, t, ...
# This is the case for getEncounter, getResourceMort, getPredRate,
# getEGrowth, getERepro, getEReproAndGrowth, getRDI
# Some have these arguments and in addition the effort argument:
#   getRates, getMort
# Some functions also accept a MizerSim and return the rate at the saved times:
# getEncounter(object, n, n_pp, n_other, time_range, drop, ...)
# getFeedingLevel(object, n, n_pp, n_other, time_range, drop, ...)
# getEReproAndGrowth(object, n, n_pp, n_other, time_range, drop, ...)
# getPredRate(object, n, n_pp, n_other, time_range, drop, ...)
# getPredMort(object, n, n_pp, n_other, time_range, drop, ...)
# getFMort(object, effort, time_range, drop, n, n_pp, n_other, t, ...)
# getFMortGear(object, effort, time_range, n, n_pp, n_other, t, ...)
# getMort(object, n, n_pp, n_other, effort, time_range, drop, ...)
# getERepro(object, n, n_pp, n_other, time_range, drop, ...)
# getEGrowth(object, n, n_pp, n_other, time_range, drop, ...)
# getDiffusion(object, n, n_pp, n_other, t, time_range, drop, ...)
# getRDI(object, n, n_pp, n_other, time_range, ...)
# getRDD(object, n, n_pp, n_other, time_range, ...)
# getFlux(object, n, n_pp, n_other, time_range, drop, ...)

#' Get all rates
#'
#' Calls other rate functions in sequence and collects the results in a list.
#' The rates returned are encounter, feeding level, energy for growth and
#' reproduction, predation rate, predation mortality, and resource mortality.
#' The purpose of this function is to provide a convenient way to get all the
#' rates at once, and to ensure that they are all calculated at the same time
#' step with the same inputs. The rates are returned in a list with the same
#' names as the rate functions that calculate them, so for example the encounter
#' rate is returned in the list element named "encounter" and is calculated with
#' the getEncounter() function.
#'
#' When mizer needs to calculate the rates during a simulation it does not use
#' this function but instead the faster `projectRates()`.
#'
#' @inheritParams mizerRates
#' @export
#' @family rate functions
#' @examples
#' rates <- getRates(NS_params)
#' names(rates)
#' identical(rates$encounter, getEncounter(NS_params))
getRates <- function(params, n = initialN(params),
                     n_pp = initialNResource(params),
                     n_other = initialNOther(params),
                     effort, t = 0, ...) {
    UseMethod("getRates")
}
#' @export
getRates.MizerParams <- function(params, n = initialN(params),
                     n_pp = initialNResource(params),
                     n_other = initialNOther(params),
                     effort, t = 0, ...) {
    params <- validParams(params)
    if (missing(effort)) {
        effort <- params@initial_effort
    }

    rates_fns <- projectRateFunctions(params)
    rates_fns$Rates(
        params, n = n, n_pp = n_pp, n_other = n_other,
        t = t, effort = effort,
        rates_fns = rates_fns, ...)
}

#' Get encounter rate
#'
#' Returns the rate at which a predator of species \eqn{i} and
#' weight \eqn{w} encounters food (grams/year).
#'
#' @inherit mizerEncounter details sections
#'
#' @template param_object_dots
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (predator species x predator
#'   size) with the encounter rates.
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x predator
#'   species x predator size) with the encounter rates at every time step.
#'   If `drop = TRUE` then dimensions of length 1 will be removed.
#' @export
#' @family rate functions
#' @examples
#' encounter <- getEncounter(NS_params)
#' str(encounter)
getEncounter <- function(object, ...) {
    UseMethod("getEncounter")
}
#' @export
getEncounter.MizerParams <- function(object, n = initialN(object),
                         n_pp = initialNResource(object),
                         n_other = initialNOther(object),
                         t = 0, ...) {
    params <- validParams(object)
    assert_that(is.array(n),
                is.numeric(n_pp),
                is.list(n_other),
                is.number(t),
                identical(dim(n), dim(params@initial_n)),
                identical(length(n_pp), length(params@initial_n_pp)),
                identical(length(n_other), length(params@initial_n_other))
    )
    if (usesExtensionDispatch(params)) {
        encounter <- projectEncounter(params, n = n, n_pp = n_pp,
                                      n_other = n_other, t = t, ...)
    } else {
        encounter_fn <- get(params@rates_funcs$Encounter)
        encounter <- encounter_fn(params, n = n, n_pp = n_pp,
                                  n_other = n_other, t = t, ...)
    }
    ArraySpeciesBySize(encounter, value_name = "Encounter rate",
             units = "g/year", params = params)
}

#' @export
getEncounter.MizerSim <- function(object, n, n_pp, n_other, t, ...,
                                  time_range, drop = FALSE) {
    sim <- object
    if (missing(time_range) && !missing(t)) time_range <- t
    sim_size_rate(sim, time_range, drop, target = "Encounter",
                  slot = "encounter", value_name = "Encounter rate",
                  units = "g/year", ...)
}


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
#'
#' @return A time x species x size array, possibly with dimensions dropped.
#' @keywords internal
get_species_size_rate_from_sim <- function(sim, time_range, drop,
                                           rate_fun, value_name,
                                           units = NULL) {
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
                                          params = sim@params)
    } else if (is.matrix(result) &&
               names(dimnames(result))[[1]] == "sp") {
        result <- ArraySpeciesBySize(result, value_name = value_name,
                                     units = units, params = sim@params)
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
                          value_name, units = NULL,
                          use_sim_effort = FALSE, ...) {
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
        value_name = value_name, units = units)
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


#' Get feeding level
#'
#' Returns the feeding level.
#' By default this function uses [mizerFeedingLevel()] to calculate
#' the feeding level, but this can be overruled via [setRateFunction()].
#'
#' @inherit mizerFeedingLevel details sections
#'
#' @template param_object_dots
#'
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (predator species x predator
#'   size) with the feeding level.
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x predator
#'   species x predator size) with the feeding level at every time step.
#'   If `drop = TRUE` then dimensions of length 1 will be removed.
#'
#' @export
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Get initial feeding level
#' fl <- getFeedingLevel(params)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at all saved time steps
#' fl <- getFeedingLevel(sim)
#' # Get the feeding level for years 15 - 20
#' fl <- getFeedingLevel(sim, time_range = c(15, 20))
#' }
getFeedingLevel <- function(object, ...) {
    UseMethod("getFeedingLevel")
}

#' @export
getFeedingLevel.MizerParams <- function(object, n, n_pp, n_other,
                            time_range, drop = FALSE, ...) {
    params <- validParams(object)
    if (missing(time_range)) time_range <- 0
    t <- min(time_range)
    if (missing(n)) n <- params@initial_n
    if (missing(n_pp)) n_pp <- params@initial_n_pp
    if (missing(n_other)) n_other <- params@initial_n_other
    # calculate feeding level
    encounter <- getEncounter(params, n, n_pp, n_other, t = t)
    if (usesExtensionDispatch(params)) {
        feeding_level <- projectFeedingLevel(
            params, n = n, n_pp = n_pp, n_other = n_other,
            encounter = encounter, t = t)
    } else {
        f <- get(params@rates_funcs$FeedingLevel)
        feeding_level <- f(params, n = n, n_pp = n_pp, n_other = n_other,
                           encounter = encounter, t = t)
    }
    return(ArraySpeciesBySize(feeding_level, value_name = "Feeding level",
                     params = params))
}

#' @export
getFeedingLevel.MizerSim <- function(object, n, n_pp, n_other,
                            time_range, drop = FALSE, ...) {
    sim <- object
    sim_size_rate(sim, time_range, drop, target = "FeedingLevel",
                  slot = "feeding_level", value_name = "Feeding level", ...)
}


#' Get critical feeding level
#'
#' The critical feeding level is the feeding level at which the food intake is
#' just high enough to cover the metabolic costs, with nothing left over for
#' growth or reproduction.
#'
#' @param params A MizerParams object
#' @return An `ArraySpeciesBySize` object (species x size) with the critical feeding level
#' @export
#' @examples
#' \donttest{
#' str(getFeedingLevel(NS_params))
#' }
getCriticalFeedingLevel <- function(params) {
    UseMethod("getCriticalFeedingLevel")
}
#' @export
getCriticalFeedingLevel.MizerParams <- function(params) {
    params <- validParams(params)
    result <- params@metab / params@intake_max / params@species_params$alpha
    ArraySpeciesBySize(result, value_name = "Critical feeding level",
             params = params)
}


#' Get energy rate available for reproduction and growth
#'
#' Calculates the energy rate \eqn{E_{r.i}(w)} (grams/year) available for
#' reproduction and growth after metabolism and movement have been accounted
#' for.
#'
#' @inherit mizerEReproAndGrowth details sections
#'
#' @template param_object_dots
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (species x size) with the
#'   energy rate \eqn{E_{r.i}(w)} available for growth and reproduction
#'   (grams/year).
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species x
#'   size) with the energy rate at every time step. If `drop = TRUE` then
#'   dimensions of length 1 will be removed.
#' @export
#' @seealso The part of this energy rate that is invested into growth is
#'   calculated with [getEGrowth()] and the part that is invested into
#'   reproduction is calculated with [getERepro()].
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' e <- getEReproAndGrowth(params, n = N(sim)[15, , ],
#'                         n_pp = NResource(sim)[15, ], t = 15)
#' # Rate at this time for Sprat of size 2g
#' e["Sprat", "2"]
#' }
getEReproAndGrowth <- function(object, ...) {
    UseMethod("getEReproAndGrowth")
}

#' @export
getEReproAndGrowth.MizerParams <- function(object, n = initialN(object),
                               n_pp = initialNResource(object),
                               n_other = initialNOther(object),
                               t = 0, ...) {
    params <- validParams(object)
    encounter <- getEncounter(params, n, n_pp, n_other, t = t)
    feeding_level <- getFeedingLevel(params, n, n_pp, n_other, time_range = t)
    if (usesExtensionDispatch(params)) {
        e <- projectEReproAndGrowth(
            params, n = n, n_pp = n_pp, n_other = n_other, t = t,
            encounter = encounter, feeding_level = feeding_level)
    } else {
        f <- get(params@rates_funcs$EReproAndGrowth)
        e <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
               encounter = encounter, feeding_level = feeding_level)
    }
    ArraySpeciesBySize(e, value_name = "Energy for growth and reproduction",
             units = "g/year", params = params)
}

#' @export
getEReproAndGrowth.MizerSim <- function(object, n, n_pp, n_other, t, ...,
                                        time_range, drop = FALSE) {
    sim <- object
    if (missing(time_range) && !missing(t)) time_range <- t
    sim_size_rate(sim, time_range, drop, target = "EReproAndGrowth",
                  slot = "e",
                  value_name = "Energy for growth and reproduction",
                  units = "g/year", ...)
}

#' Get predation rate
#'
#' Calculates the potential rate (in units 1/year) at which a prey individual of
#' a given size \eqn{w} is killed by predators from species \eqn{j}. In formulas
#' \deqn{{\tt pred\_rate}_j(w_p) = \int \phi_j(w,w_p) (1-f_j(w))
#'   \gamma_j(w) N_j(w) \, dw.}{pred_rate_j(w_p) = \int\phi_i(w,w_p) (1-f_i(w))
#'   \gamma_i(w) N_i(w) dw.}
#' This potential rate is used in [getPredMort()] to
#' calculate the realised predation mortality rate on the prey individual.
#'
#' @inherit mizerPredRate details sections
#'
#' @template param_object_dots
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (predator species x prey
#'   size), where the prey size runs over fish community plus resource spectrum.
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x predator
#'   species x prey size) with the predation rates at every time step.
#'   If `drop = TRUE` then dimensions of length 1 will be removed.
#' @export
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Predation rate in initial state
#' pred_rate <- getPredRate(params)
#' str(pred_rate)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' pred_rate <- getPredRate(params, n = N(sim)[15, , ],
#'                          n_pp = NResource(sim)[15, ], t = 15)
#' }
getPredRate <- function(object, ...) {
    UseMethod("getPredRate")
}

#' @export
getPredRate.MizerParams <- function(object, n = initialN(object),
                        n_pp = initialNResource(object),
                        n_other = initialNOther(object),
                        t = 0, ...) {
    params <- validParams(object)
    feeding_level <- getFeedingLevel(params, n = n, n_pp = n_pp,
                                     n_other = n_other, time_range = t)
    if (usesExtensionDispatch(params)) {
        pred_rate <- projectPredRate(
            params, n = n, n_pp = n_pp, n_other = n_other, t = t,
            feeding_level = feeding_level)
    } else {
        f <- get(params@rates_funcs$PredRate)
        pred_rate <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                       feeding_level = feeding_level)
    }
    dimnames(pred_rate) <- list(sp = dimnames(params@initial_n)$sp,
                                w_prey = as.character(signif(params@w_full, 3)))
    ArraySpeciesBySize(pred_rate, value_name = "Predation rate",
                       units = "1/year", params = params)
}

#' @export
getPredRate.MizerSim <- function(object, n, n_pp, n_other, t, ...,
                                 time_range, drop = FALSE) {
    sim <- object
    if (missing(time_range) && !missing(t)) time_range <- t
    sim_size_rate(sim, time_range, drop, target = "PredRate",
                  slot = "pred_rate", value_name = "Predation rate",
                  units = "1/year", ...)
}

#' Get total predation mortality rate
#'
#' Calculates the total predation mortality rate \eqn{\mu_{p,i}(w_p)} (in units
#' of 1/year) on each prey species by prey size:
#' \deqn{\mu_{p.i}(w_p) = \sum_j {\tt pred\_rate}_j(w_p)\, \theta_{ji}.}{
#'   \mu_{p.i}(w_p) = \sum_j pred_rate_j(w_p) \theta_{ji}.}
#' The predation rate `pred_rate` is returned by [getPredRate()].
#'
#' @inherit mizerPredMort details sections
#' @template param_object_dots
#'
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (prey species x prey size)
#'   with the predation mortality rates.
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x prey species
#'   x prey size) with the predation mortality at every time step.
#'   If `drop = TRUE` then dimensions of length 1 will be removed.
#' @family rate functions
#' @export
#' @examples
#' \donttest{
#' params <- NS_params
#' # Predation mortality in initial state
#' M2 <- getPredMort(params)
#' str(M2)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get predation mortality at one time step
#' M2 <- getPredMort(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ])
#' # Get predation mortality at all saved time steps
#' M2 <- getPredMort(sim)
#' str(M2)
#' # Get predation mortality over the years 15 - 20
#' M2 <- getPredMort(sim, time_range = c(15, 20))
#' }
getPredMort <- function(object, ...) {
    UseMethod("getPredMort")
}

#' @export
getPredMort.MizerParams <- function(object, n, n_pp, n_other,
                        time_range, drop = TRUE, ...) {
    params <- validParams(object)
    if (missing(n)) n <- params@initial_n
    if (missing(n_pp)) n_pp <- params@initial_n_pp
    if (missing(n_other)) n_other <- params@initial_n_other
    if (missing(time_range)) time_range <- 0
    t <- min(time_range)

    pred_rate <- getPredRate(params, n = n, n_pp = n_pp, n_other = n_other,
                             t = t)
    if (usesExtensionDispatch(params)) {
        pred_mort <- projectPredMort(
            params, n = n, n_pp = n_pp, n_other = n_other, t = t,
            pred_rate = pred_rate)
    } else {
        f <- get(params@rates_funcs$PredMort)
        pred_mort <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                       pred_rate = pred_rate)
    }
    ArraySpeciesBySize(pred_mort, value_name = "Predation mortality",
             units = "1/year", params = params)
}

#' @export
getPredMort.MizerSim <- function(object, n, n_pp, n_other,
                        time_range, drop = TRUE, ...) {
    sim <- object
    sim_size_rate(sim, time_range, drop, target = "PredMort",
                  slot = "pred_mort", value_name = "Predation mortality",
                  units = "1/year", ...)
}

#' Alias for `getPredMort()`
#'
#' `r lifecycle::badge("deprecated")`
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getPredMort
#' @export
#' @concept deprecated
getM2 <- getPredMort


#' Get predation mortality rate for resource
#'
#' Calculates the predation mortality rate \eqn{\mu_p(w)} on the resource
#' spectrum by resource size (in units 1/year).
#'
#' @inherit mizerResourceMort
#' @inheritParams mizerRates
#'
#' @return A vector of mortality rate by resource size.
#' @family rate functions
#' @export
#' @examples
#' \donttest{
#' params <- NS_params
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get resource mortality at one time step
#' getResourceMort(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ])
#' }
getResourceMort <-
    function(params, n = initialN(params),
             n_pp = initialNResource(params),
             n_other = initialNOther(params),
             t = 0, ...) {
    UseMethod("getResourceMort")
}

#' @export
getResourceMort.MizerParams <- function(params, n = initialN(params),
             n_pp = initialNResource(params),
             n_other = initialNOther(params),
             t = 0, ...) {
    params <- validParams(params)

    pred_rate <- getPredRate(params, n = n, n_pp = n_pp, n_other = n_other,
                             t = t)
    if (usesExtensionDispatch(params)) {
        mort <- projectResourceMort(
            params, n = n, n_pp = n_pp, n_other = n_other,
            t = t, pred_rate = pred_rate)
    } else {
        f <- get(params@rates_funcs$ResourceMort)
        mort <- f(params, n = n, n_pp = n_pp, n_other = n_other,
                  t = t, pred_rate = pred_rate)
    }
    names(mort) <- names(params@initial_n_pp)
    ArrayResourceBySize(mort, value_name = "Resource mortality",
                        units = "1/year", params = params)
}

#' Alias for `getResourceMort()`
#'
#' `r lifecycle::badge("deprecated")`
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getResourceMort
#' @export
#' @concept deprecated
getM2Background <- getResourceMort


#' Get the fishing mortality by time, gear, species and size
#'
#' Calculates the fishing mortality rate \eqn{F_{g,i,w}} by gear, species and
#' size and possibly time (in units 1/year).
#'
#' @param object A \linkS4class{MizerParams} or \linkS4class{MizerSim} object.
#' @param ... Additional arguments that depend on the class of `object`.
#'
#'   **For a \linkS4class{MizerParams} object:**
#'   \describe{
#'     \item{`effort`}{The effort for each fishing gear. See notes below.
#'       Defaults to the initial effort stored in `object`.}
#'     \item{`n`}{A matrix of species abundances (species x size). Defaults to
#'       the initial abundances stored in `object`.}
#'     \item{`n_pp`}{A vector of the resource abundance by size. Defaults to the
#'       initial resource abundance stored in `object`.}
#'     \item{`n_other`}{A named list of the abundances of other dynamical
#'       components. Defaults to the initial values stored in `object`.}
#'     \item{`t`}{The time for which to do the calculation. Defaults to 0.}
#'   }
#'
#'   **For a \linkS4class{MizerSim} object:**
#'   \describe{
#'     \item{`time_range`}{Subset the returned fishing mortalities by time. The
#'       time range is either a vector of values, a vector of min and max time,
#'       or a single value. Defaults to the whole time range of the
#'       simulation.}
#'   }
#'
#' @return An array. If the effort argument has a time dimension, or a
#'   `MizerSim` is passed in, the output array has four dimensions (time x
#'   gear x species x size). If the effort argument does not have a time
#'   dimension (i.e. it is a vector or a single numeric), the output array has
#'   three dimensions (gear x species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#'
#' The `effort` argument is only used if a `MizerParams` object is
#' passed in. The `effort` argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' `effort` argument must be the same the same as in the `MizerParams`
#' object. If the `effort` argument is not supplied, its value is taken
#' from the `@initial_effort` slot in the params object.
#'
#' If the object argument is of class `MizerSim` then the effort slot of
#' the `MizerSim` object is used and the `effort` argument is not
#' used.
#' @export
#' @family rate functions
#' @examples
#' \donttest{
#' params <-NS_params
#' # Get the fishing mortality in initial state
#' F <- getFMortGear(params, effort = 1)
#' str(F)
#' # Get the initial fishing mortality when effort is different
#' # between the four gears:
#' F <- getFMortGear(params, effort = c(0.5, 1, 1.5, 0.75))
#' # Get the fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20, 4))
#' effort[, 1] <- seq(from=0, to = 1, length = 20)
#' effort[, 2] <- seq(from=1, to = 0.5, length = 20)
#' effort[, 3] <- seq(from=1, to = 2, length = 20)
#' effort[, 4] <- seq(from=2, to = 1, length = 20)
#' F <- getFMortGear(params, effort = effort)
#' str(F)
#' # Get the fishing mortality using the effort already held in a MizerSim object.
#' sim <- project(params, t_max = 20, effort = 0.5)
#' F <- getFMortGear(sim)
#' F <- getFMortGear(sim, time_range = c(10, 20))
#' }
#'
getFMortGear <- function(object, ...) {
    UseMethod("getFMortGear")
}

#' @export
getFMortGear.MizerParams <- function(object, effort, time_range,
                                     n, n_pp, n_other, t = 0, ...) {
    params <- validParams(object)
    if (missing(effort)) {
        effort <- params@initial_effort
    }
    if (is(effort, "numeric")) {
            no_gear <- dim(params@catchability)[1]
            # If a single value, just repeat it for all gears
            if (length(effort) == 1) {
                effort <- rep(effort, no_gear)
            }
            if (length(effort) != no_gear) {
                stop("Effort must be a single value or a vector as long as the number of gears\n")
            }

            f <- mizerFMortGear(params, effort = effort)
            dimnames(f) <- dimnames(params@selectivity)
            return(f)
        } else {
            # assuming effort is a matrix, and object is of MizerParams class
            no_gear <- dim(params@catchability)[1]
            if (dim(effort)[2] != no_gear)
                stop("Effort array must have a single value or a vector as long as the number of gears for each time step\n")
            # Make the output array - note that we put time as last dimension
            # and then aperm before returning. This is because of the order of
            # the values when we call the other getFMortGear function.
            # Fill it up by calling the other function and passing in each line
            # of the effort matrix
            out <- array(NA, dim = c(dim(params@selectivity), dim(effort)[1]),
                         dimnames = c(dimnames(params@selectivity),
                                      list(time = dimnames(effort)[[1]])))
            out[] <- apply(effort, 1, function(x) mizerFMortGear(params, x))
            out <- aperm(out, c(4, 1, 2, 3))
            return(out)
        }
}

#' @export
getFMortGear.MizerSim <- function(object, effort, time_range,
                                  n, n_pp, n_other, t = 0, ...) {
    sim <- object
    time_elements <- get_sim_rate_time_elements(sim, time_range)
    # Validate the params once rather than at every saved time step.
    params <- validParams(sim@params)

    # Work slice by slice, like the other MizerSim rate methods, so any future
    # state-aware gear mortality calculation receives the matching simulation
    # state and time rather than only the effort matrix.
    f_mort_gear <- plyr::aaply(which(time_elements), 1, function(time_idx) {
        slice <- get_sim_rate_slice(sim, time_idx)
        getFMortGear(params, effort = slice$effort, n = slice$n,
                     n_pp = slice$n_pp, n_other = slice$n_other,
                     t = slice$t, ...)
    }, .drop = FALSE)
    names(dimnames(f_mort_gear))[[1]] <- "time"
    f_mort_gear
}

#' Get the total fishing mortality rate from all fishing gears by time, species
#' and size.
#'
#' Calculates the total fishing mortality  (in units 1/year) from all gears by
#' species and size and possibly time. See `setFishing()` for details of
#' how fishing gears are set up.
#'
#' The total fishing mortality is just the sum of the fishing mortalities
#' imposed by each gear, \eqn{F_i(w)=\sum_g F_{g,i,w}}.
#' The fishing mortality for each gear is obtained as catchability x
#' selectivity x effort.
#'
#' @param object A \linkS4class{MizerParams} or \linkS4class{MizerSim} object.
#' @param ... Additional arguments that depend on the class of `object`.
#'
#'   **For a \linkS4class{MizerParams} object:**
#'   \describe{
#'     \item{`effort`}{The effort of each fishing gear. See notes below.
#'       Defaults to the initial effort stored in `object`.}
#'     \item{`n`}{A matrix of species abundances (species x size). Defaults to
#'       the initial abundances stored in `object`.}
#'     \item{`n_pp`}{A vector of the resource abundance by size. Defaults to the
#'       initial resource abundance stored in `object`.}
#'     \item{`n_other`}{A named list of the abundances of other dynamical
#'       components. Defaults to the initial values stored in `object`.}
#'     \item{`t`}{The time for which to do the calculation. Defaults to 0.}
#'   }
#'
#'   **For a \linkS4class{MizerSim} object:**
#'   \describe{
#'     \item{`time_range`}{Subset the returned fishing mortalities by time. The
#'       time range is either a vector of values, a vector of min and max time,
#'       or a single value. Defaults to the whole time range of the
#'       simulation.}
#'     \item{`drop`}{Should dimensions of length 1 be dropped, e.g. if your
#'       community only has one species it might make presentation of results
#'       easier. Defaults to `TRUE`.}
#'   }
#'
#' @return
#' * `MizerParams` with vector effort: An `ArraySpeciesBySize` object
#'   (species x size) with the fishing mortality rates.
#' * `MizerParams` with time-dimensioned effort or `MizerSim`: An
#'   `ArrayTimeBySpeciesBySize` object (time x species x size).
#'
#' The `effort` argument is only used if a `MizerParams` object is
#' passed in. The `effort` argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' `effort` argument must be the same as in the `MizerParams`
#' object.
#'
#' If the object argument is of class `MizerSim` then the effort slot of
#' the `MizerSim` object is used and the `effort` argument is not
#' used.
#'
#' @inheritSection mizerFMort Your own fishing mortality function
#'
#' @export
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Get the total fishing mortality in the initial state
#' F <- getFMort(params, effort = 1)
#' str(F)
#' # Get the initial total fishing mortality when effort is different
#' # between the four gears:
#' F <- getFMort(params, effort = c(0.5,1,1.5,0.75))
#' # Get the total fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20,4))
#' effort[, 1] <- seq(from = 0, to = 1, length = 20)
#' effort[, 2] <- seq(from = 1, to = 0.5, length = 20)
#' effort[, 3] <- seq(from = 1, to = 2, length = 20)
#' effort[, 4] <- seq(from = 2, to = 1, length = 20)
#' F <- getFMort(params, effort = effort)
#' str(F)
#' # Get the total fishing mortality using the effort already held in a
#' # MizerSim object.
#' sim <- project(params, t_max = 20, effort = 0.5)
#' F <- getFMort(sim)
#' F <- getFMort(sim, time_range = c(10, 20))
#' }
getFMort <- function(object, ...) {
    UseMethod("getFMort")
}

#' @export
getFMort.MizerParams <- function(object, effort, time_range, drop = TRUE,
                                 n, n_pp, n_other, t, ...) {
    params <- validParams(object)
    if (missing(effort)) {
        effort <- params@initial_effort
    }
    if (missing(t)) {
        # Legacy support: before `getFMort()` accepted `t` directly,
        # `time_range` was used to set the time passed to custom FMort
        # functions. New code should pass `t`.
        if (missing(time_range)) time_range <- 0
        t <- min(time_range)
    }
    if (missing(n)) n <- params@initial_n
    if (missing(n_pp)) n_pp <- params@initial_n_pp
    if (missing(n_other)) n_other <- params@initial_n_other
    no_gears <- dim(params@catchability)[[1]]
    f <- if (usesExtensionDispatch(params)) {
        projectFMort
    } else {
        get(params@rates_funcs$FMort)
    }
    if (length(dim(effort)) == 2) {
        times <- dimnames(effort)$time
        f_mort <- array(0,
                        dim = c(dim(effort)[[1]], dim(params@initial_n)),
                            dimnames = c(list(time = times),
                                         dimnames(params@initial_n)))
        times <- as.numeric(times)
        for (i in seq_len(dim(effort)[1])) {
            args <- list(
                params = params, n = n, n_pp = n_pp, n_other = n_other,
                effort = effort[i, ], t = times[i],
                e_growth = getEGrowth(params, n = n, n_pp = n_pp,
                                      n_other = n_other, t = times[i]),
                pred_mort = getPredMort(params, n = n, n_pp = n_pp,
                                        n_other = n_other,
                                        time_range = times[i]))
            f_mort[i, , ] <- do.call(f, c(args, list(...)))
        }
        return(f_mort)
    } else if (length(effort) <= 1) {
        args <- list(
            params = params, n = n, n_pp = n_pp, n_other = n_other,
            effort = rep(effort, no_gears), t = t,
            e_growth = getEGrowth(params, n = n, n_pp = n_pp,
                                  n_other = n_other, t = t),
            pred_mort = getPredMort(params, n = n, n_pp = n_pp,
                                    n_other = n_other, time_range = t))
        fmort <- do.call(f, c(args, list(...)))
        fmort <- ArraySpeciesBySize(fmort, value_name = "Fishing mortality",
                           units = "1/year", params = params)
        return(fmort)
    } else if (length(effort) == no_gears) {
        args <- list(
            params = params, n = n, n_pp = n_pp, n_other = n_other,
            effort = effort, t = t,
            e_growth = getEGrowth(params, n = n, n_pp = n_pp,
                                  n_other = n_other, t = t),
            pred_mort = getPredMort(params, n = n, n_pp = n_pp,
                                    n_other = n_other, time_range = t))
        fmort <- do.call(f, c(args, list(...)))
        fmort <- ArraySpeciesBySize(fmort, value_name = "Fishing mortality",
                           units = "1/year", params = params)
        return(fmort)
    } else {
        stop("Invalid effort argument")
    }
}
#'
#' @export
getFMort.MizerSim <- function(object, effort, time_range, drop = TRUE,
                              n, n_pp, n_other, t, ...) {
    sim <- object
    sim_size_rate(sim, time_range, drop, target = "FMort", slot = "f_mort",
                  value_name = "Fishing mortality", units = "1/year",
                  use_sim_effort = TRUE, ...)
}


#' Get total mortality rate
#'
#' Calculates the total mortality rate \eqn{\mu_i(w)}  (in units 1/year) on each
#' species by size from predation mortality, background mortality and fishing
#' mortality for a single time step.
#'
#' @inherit mizerMort details sections
#' @param object A \linkS4class{MizerParams} or \linkS4class{MizerSim} object.
#' @param ... Additional arguments that depend on the class of `object`.
#'
#'   **For a \linkS4class{MizerParams} object:**
#'   \describe{
#'     \item{`n`}{A matrix of species abundances (species x size). Defaults to
#'       the initial abundances stored in `object`.}
#'     \item{`n_pp`}{A vector of the resource abundance by size. Defaults to the
#'       initial resource abundance stored in `object`.}
#'     \item{`n_other`}{A named list of the abundances of other dynamical
#'       components. Defaults to the initial values stored in `object`.}
#'     \item{`effort`}{A numeric vector of the effort by gear or a single
#'       numeric effort value which is used for all gears. Defaults to the
#'       initial effort stored in `object`.}
#'     \item{`t`}{The time for which to do the calculation. Defaults to 0.}
#'   }
#'
#'   **For a \linkS4class{MizerSim} object:**
#'   \describe{
#'     \item{`time_range`}{The time range over which to return the rates. Either
#'       a vector of values, a vector of min and max time, or a single value.
#'       Defaults to the whole time range of the simulation.}
#'     \item{`drop`}{If `TRUE` then any dimension of length 1 is removed from the
#'       returned array.}
#'   }
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (species x size) with the
#'   total mortality rates.
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species x
#'   size) with the total mortality rates at every time step. If `drop = TRUE`
#'   then dimensions of length 1 will be removed.
#'
#' @export
#' @seealso [getPredMort()], [getFMort()]
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the total mortality at a particular time step
#' mort <- getMort(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ],
#'                 t = 15, effort = 0.5)
#' # Mortality rate at this time for Sprat of size 2g
#' mort["Sprat", "2"]
#' }
getMort <- function(object, ...) {
    UseMethod("getMort")
}

#' @export
getMort.MizerParams <- function(object,
                    n = initialN(object),
                    n_pp = initialNResource(object),
                    n_other = initialNOther(object),
                    effort = getInitialEffort(object),
                    t = 0, ...) {
    params <- validParams(object)
    rates_fns <- projectRateFunctions(params)
    r <- mizer_rates_subset(params, n = n, n_pp = n_pp, n_other = n_other,
                             t = t, effort = effort,
                             rates_fns = rates_fns, targets = "Mort", ...)
    return(ArraySpeciesBySize(r$mort, value_name = "Total mortality",
                     units = "1/year", params = params))
}

#' @export
getMort.MizerSim <- function(object, n, n_pp, n_other, effort, t, ...,
                             time_range, drop = TRUE) {
    sim <- object
    if (missing(time_range) && !missing(t)) time_range <- t
    sim_size_rate(sim, time_range, drop, target = "Mort", slot = "mort",
                  value_name = "Total mortality", units = "1/year",
                  use_sim_effort = TRUE, ...)
}

#' Alias for `getMort()`
#'
#' `r lifecycle::badge("deprecated")`
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getMort
#' @export
#' @concept deprecated
getZ <- getMort


#' Get energy rate available for reproduction
#'
#' Calculates the energy rate (grams/year) available for reproduction after
#' growth and metabolism have been accounted for.
#'
#' @inherit mizerERepro details sections
#' @template param_object_dots
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (species x size) holding
#'   \deqn{\psi_i(w)\max(0, E_{r.i}(w))}
#'   where \eqn{E_{r.i}(w)} is the rate at which energy becomes available for
#'   growth and reproduction, calculated with [getEReproAndGrowth()],
#'   and \eqn{\psi_i(w)} is the proportion of this energy that is used for
#'   reproduction. Negative values of \eqn{E_{r.i}(w)} are clipped to 0 before
#'   multiplying by \eqn{\psi_i(w)}. This proportion is taken from the `params`
#'   object and is set with [setReproduction()].
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species x
#'   size) with the energy for reproduction at every time step. If `drop = TRUE`
#'   then dimensions of length 1 will be removed.
#' @export
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the rate at a particular time step
#' erepro <- getERepro(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
#' # Rate at this time for Sprat of size 2g
#' erepro["Sprat", "2"]
#' }
getERepro <- function(object, ...) {
    UseMethod("getERepro")
}

#' @export
getERepro.MizerParams <- function(object, n = initialN(object),
                      n_pp = initialNResource(object),
                      n_other = initialNOther(object),
                      t = 0, ...) {
    params <- validParams(object)
    e <- getEReproAndGrowth(params, n = n, n_pp = n_pp, n_other = n_other,
                            t = t)
    if (usesExtensionDispatch(params)) {
        erepro <- projectERepro(params, n = n, n_pp = n_pp, n_other = n_other,
                                t = t, e = e)
    } else {
        f <- get(params@rates_funcs$ERepro)
        erepro <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                    e = e)
    }
    ArraySpeciesBySize(erepro, value_name = "Energy for reproduction",
             units = "g/year", params = params)
}

#' @export
getERepro.MizerSim <- function(object, n, n_pp, n_other, t, ...,
                               time_range, drop = FALSE) {
    sim <- object
    if (missing(time_range) && !missing(t)) time_range <- t
    sim_size_rate(sim, time_range, drop, target = "ERepro", slot = "e_repro",
                  value_name = "Energy for reproduction", units = "g/year", ...)
}

#' Alias for `getERepro()`
#'
#' `r lifecycle::badge("deprecated")`
#' An alias provided for backward compatibility with mizer version <= 1.0
#' @inherit getERepro
#' @export
#' @concept deprecated
getESpawning <- getERepro


#' Get energy rate available for growth
#'
#' Calculates the energy rate \eqn{g_i(w)} (grams/year) available by species and
#' size for growth after metabolism, movement and reproduction have been
#' accounted for.
#'
#' @inherit mizerEGrowth details sections
#' @template param_object_dots
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (species x size) with the
#'   somatic growth rates (grams/year).
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species x
#'   size) with the growth rates at every time step. If `drop = TRUE` then
#'   dimensions of length 1 will be removed.
#' @export
#' @seealso [getERepro()], [getEReproAndGrowth()]
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' growth <- getEGrowth(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
#' # Growth rate at this time for Sprat of size 2g
#' growth["Sprat", "2"]
#' }
getEGrowth <- function(object, ...) {
    UseMethod("getEGrowth")
}

#' @export
getEGrowth.MizerParams <- function(object, n = initialN(object),
                       n_pp = initialNResource(object),
                       n_other = initialNOther(object),
                       t = 0, ...) {
    params <- validParams(object)
    e_repro <- getERepro(params, n = n, n_pp = n_pp, n_other = n_other,
                         t = t)
    e <- getEReproAndGrowth(params, n = n, n_pp = n_pp,
                            n_other = n_other, t = t)
    if (usesExtensionDispatch(params)) {
        g <- projectEGrowth(params, n = n, n_pp = n_pp, n_other = n_other,
                            t = t, e_repro = e_repro, e = e)
    } else {
        f <- get(params@rates_funcs$EGrowth)
        g <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
               e_repro = e_repro, e = e)
    }
    ArraySpeciesBySize(g, value_name = "Growth rate",
             units = "g/year", params = params)
}

#' @export
getEGrowth.MizerSim <- function(object, n, n_pp, n_other, t, ...,
                                time_range, drop = FALSE) {
    sim <- object
    if (missing(time_range) && !missing(t)) time_range <- t
    sim_size_rate(sim, time_range, drop, target = "EGrowth", slot = "e_growth",
                  value_name = "Growth rate", units = "g/year", ...)
}


#' Get density independent rate of egg production
#'
#' Calculates the density-independent rate of total egg production
#' \eqn{R_{di}}{R_di} (units 1/year) before density dependence, by species.
#'
#' @inherit mizerRDI details sections
#' @template param_object_dots_nodrop
#'
#' @return
#' * `MizerParams`: A numeric vector the length of the number of species.
#' * `MizerSim`: An `ArrayTimeBySpecies` object (time x species).
#' @export
#' @seealso [getRDD()]
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the density-independent reproduction rate at a particular time step
#' getRDI(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
#' }
getRDI <- function(object, ...) {
    UseMethod("getRDI")
}

#' @export
getRDI.MizerParams <- function(object, n = initialN(object),
                   n_pp = initialNResource(object),
                   n_other = initialNOther(object),
                   t = 0, ...) {
    params <- validParams(object)
    e_repro <- getERepro(params, n = n, n_pp = n_pp, n_other = n_other,
                         t = t)
    e_growth <- getEGrowth(params, n = n, n_pp = n_pp, n_other = n_other,
                           t = t)
    diffusion <- getDiffusion(params, n = n, n_pp = n_pp, n_other = n_other,
                              t = t)
    mort <- getMort(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
    if (usesExtensionDispatch(params)) {
        rdi <- projectRDI(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                          e_repro = e_repro, e_growth = e_growth, mort = mort,
                          diffusion = diffusion)
    } else {
        f <- get(params@rates_funcs$RDI)
        rdi <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                 e_repro = e_repro, e_growth = e_growth, mort = mort,
                 diffusion = diffusion)
    }
    names(rdi) <- params@species_params$species
    rdi
}

#' @export
getRDI.MizerSim <- function(object, n, n_pp, n_other, t = 0,
                            time_range, ...) {
    sim <- object
    sim_species_rate(sim, time_range, target = "RDI", slot = "rdi",
                     value_name = "Density-independent reproduction rate",
                     units = "1/year", use_sim_effort = TRUE, ...)
}


#' Get density dependent reproduction rate
#'
#' Calculates the density dependent rate of egg production \eqn{R_i} (units
#' 1/year) for each species. This is the flux entering the smallest size class
#' of each species. The density dependent rate is the density independent
#' rate obtained with [getRDI()] after it has been put through the
#' density dependence function. This is the Beverton-Holt function
#' [BevertonHoltRDD()] by default, but this can be changed. See
#' [setReproduction()] for more details.
#'
#' @param object A \linkS4class{MizerParams} or \linkS4class{MizerSim} object.
#' @param ... Additional arguments that depend on the class of `object`.
#'
#'   **For a \linkS4class{MizerParams} object:**
#'   \describe{
#'     \item{`n`}{A matrix of species abundances (species x size). Defaults to
#'       the initial abundances stored in `object`.}
#'     \item{`n_pp`}{A vector of the resource abundance by size. Defaults to the
#'       initial resource abundance stored in `object`.}
#'     \item{`n_other`}{A named list of the abundances of other dynamical
#'       components. Defaults to the initial values stored in `object`.}
#'     \item{`t`}{The time for which to do the calculation. Defaults to 0.}
#'     \item{`rdi`}{A vector of density-independent reproduction rates for each
#'       species. If not specified, it is calculated internally using
#'       [getRDI()].}
#'   }
#'
#'   **For a \linkS4class{MizerSim} object:**
#'   \describe{
#'     \item{`time_range`}{The time range over which to return the rates. Either
#'       a vector of values, a vector of min and max time, or a single value.
#'       Defaults to the whole time range of the simulation.}
#'   }
#'
#' @return
#' * `MizerParams`: A numeric vector the length of the number of species.
#' * `MizerSim`: An `ArrayTimeBySpecies` object (time x species).
#' @export
#' @seealso [getRDI()]
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the rate at a particular time step
#' getRDD(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
#' }
getRDD <- function(object, ...) {
    UseMethod("getRDD")
}

#' @export
getRDD.MizerParams <- function(object, n = initialN(object),
                   n_pp = initialNResource(object),
                   n_other = initialNOther(object),
                   t = 0,
                   rdi = getRDI(object, n = n, n_pp = n_pp,
                                n_other = n_other, t = t), ...) {
    params <- validParams(object)
    if (usesExtensionDispatch(params)) {
        rdd <- projectRDD(params, rdi = rdi, species_params = params@species_params,
                          t = t, ...)
    } else {
        # Avoid getting into infinite loops
        if (params@rates_funcs$RDD == "getRDD") {
            stop('"getRDD" is not a valid name for the function giving the density',
                 'dependent reproductive rate.')
        }
        f <- get(params@rates_funcs$RDD)
        rdd <- f(rdi = rdi, species_params = params@species_params,
                 params = params, t = t)
    }
    names(rdd) <- params@species_params$species
    rdd
}

#' @export
getRDD.MizerSim <- function(object, n, n_pp, n_other, t = 0,
                            rdi, time_range, ...) {
    sim <- object
    sim_species_rate(sim, time_range, target = "RDD", slot = "rdd",
                     value_name = "Density-dependent reproduction rate",
                     units = "1/year", use_sim_effort = TRUE, ...)
}

#' Get flux into size bins
#'
#' Calculates the flux \eqn{J_i(w)} (numbers/year) entering each size class
#' from the one below it. This is composed of an advective flux from somatic
#' growth and a diffusive flux from the redistribution of individuals.
#'
#' At the recruitment size, the flux is simply the recruitment rate
#' \eqn{R_{dd,i}} (see [getRDD()]). For sizes below the recruitment size
#' the flux is zero.
#'
#' The flux at weight \eqn{w} is multiplied by \eqn{w} raised to the `power`
#' given by the `power` argument, similar to the `power` argument of
#' [plotSpectra()]. The default `power = 0` returns the flux of individuals
#' (numbers/year). With `power = 1` the result is the flux of biomass
#' (grams/year).
#'
#' @template param_object_dots
#' @param power The flux at weight \eqn{w} is multiplied by \eqn{w} raised to
#'   `power`. The default `power = 0` gives the flux of individuals
#'   (numbers/year), whereas `power = 1` gives the flux of biomass
#'   (grams/year).
#' @return
#' * `MizerParams`: An `ArraySpeciesBySize` object (species x size) with the
#'   flux entering each size class. The units are `numbers/year` when
#'   `power = 0` and `g^power/year` otherwise.
#' * `MizerSim`: An `ArrayTimeBySpeciesBySize` object (time step x species x
#'   size) with the flux at every time step. If `drop = TRUE` then dimensions
#'   of length 1 will be removed.
#' @export
#' @seealso [getEGrowth()], [getRDD()]
#' @family rate functions
#' @examples
#' \donttest{
#' params <- NS_params
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project(params, t_max = 20, effort = 0.5)
#' # Get the flux at a particular time step
#' flux <- getFlux(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
#' # Flux for Sprat of size 2g
#' flux["Sprat", "2"]
#' }
getFlux <- function(object, ..., power = 0) {
    UseMethod("getFlux")
}

#' @rdname getFlux
#' @usage NULL
#' @export
getFlux.MizerParams <- function(object, n = initialN(object),
                    n_pp = initialNResource(object),
                    n_other = initialNOther(object),
                    t = 0, power = 0, ...) {
    params <- validParams(object)
    assert_that(is.number(power))

    g <- getEGrowth(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
    d <- getDiffusion(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
    rdd <- getRDD(params, n = n, n_pp = n_pp, n_other = n_other, t = t)

    flux <- flux_from_rates(params, n = n, g = g, d = d, rdd = rdd,
                            power = power)

    ArraySpeciesBySize(flux, value_name = "Flux",
             units = flux_units(power), params = params)
}

#' Assemble the flux matrix from growth, diffusion and recruitment rates
#'
#' Internal helper holding the arithmetic shared by `getFlux.MizerParams` and
#' `getFlux.MizerSim`. Keeping it separate lets the `MizerSim` method resolve the
#' rate functions once and reuse them across all saved time steps.
#'
#' @param params A valid `MizerParams` object.
#' @param n A matrix of species abundances (species x size).
#' @param g Growth rate matrix (species x size), as from [getEGrowth()].
#' @param d Diffusion rate matrix (species x size), as from [getDiffusion()].
#' @param rdd Density-dependent reproduction rate vector (one per species), as
#'   from [getRDD()].
#' @param power The flux at weight \eqn{w} is multiplied by \eqn{w} raised to
#'   `power`. The default `power = 0` leaves the flux of individuals unchanged.
#'
#' @return A plain species x size matrix of fluxes (no mizer array class).
#' @keywords internal
flux_from_rates <- function(params, n, g, d, rdd, power = 0) {
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    dw <- params@dw

    flux <- matrix(0, nrow = no_sp, ncol = no_w,
                   dimnames = dimnames(params@metab))

    idx <- 2:no_w
    idx_minus_1 <- idx - 1

    # Calculate J_{i,j} for all j > 1
    # J_{i,j} = g_{i, j-1} N_{i, j-1} - 1/2 * (d_{i, j} N_{i, j} - d_{i, j-1} N_{i, j-1}) / dw_{j-1}
    diff_term <- (d[, idx] * n[, idx] - d[, idx_minus_1] * n[, idx_minus_1]) /
        matrix(dw[idx_minus_1], nrow = no_sp, ncol = length(idx_minus_1), byrow = TRUE)

    flux[, idx] <- g[, idx_minus_1] * n[, idx_minus_1] - 0.5 * diff_term

    # Apply recruitment boundary conditions
    j_start <- params@w_min_idx
    idxs <- cbind(1:no_sp, j_start)
    flux[idxs] <- rdd

    # Zero out elements for sizes smaller than w_min_idx
    w_idx_mat <- matrix(1:no_w, nrow = no_sp, ncol = no_w, byrow = TRUE)
    mask_below <- w_idx_mat < j_start

    if (any(mask_below)) {
        flux[mask_below] <- 0
    }

    if (power != 0) {
        flux <- sweep(flux, 2, params@w^power, "*")
    }

    flux
}

#' Units string for the flux returned by [getFlux()]
#'
#' @param power The power of weight the flux was multiplied by.
#' @return A character string with the units of the flux.
#' @keywords internal
flux_units <- function(power) {
    if (power == 0) {
        "1/year"
    } else if (power == 1) {
        "g/year"
    } else {
        paste0("g^", power, "/year")
    }
}

#' @rdname getFlux
#' @usage NULL
#' @export
getFlux.MizerSim <- function(object, n, n_pp, n_other, t, power = 0, ...,
                             time_range, drop = FALSE) {
    sim <- object
    assert_that(is.number(power))
    if (missing(time_range) && !missing(t)) time_range <- t

    # Resolve the rate functions once rather than at every saved time step. The
    # flux needs only the growth, diffusion and recruitment rates (and the rates
    # those depend on), which `mizer_rates_subset()` calculates selectively.
    params <- validParams(sim@params)
    rates_fns <- projectRateFunctions(params)

    get_species_size_rate_from_sim(
        sim, time_range, drop,
        function(slice) {
            r <- mizer_rates_subset(
                params, n = slice$n, n_pp = slice$n_pp,
                n_other = slice$n_other, t = slice$t, effort = slice$effort,
                rates_fns = rates_fns,
                targets = c("EGrowth", "Diffusion", "RDD"), ...)
            flux_from_rates(params, n = slice$n,
                            g = r$e_growth, d = r$diffusion, rdd = r$rdd,
                            power = power)
        },
        value_name = "Flux", units = flux_units(power))
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
