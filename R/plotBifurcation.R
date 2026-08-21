# Bifurcation diagram over fishing effort -------------------------------------
#
# `plotBifurcation()` sweeps the fishing effort, follows the attractor of the
# full dynamics at each effort value, and plots the long-term range (envelope)
# of a chosen summary quantity. Where the steady state is stable the envelope
# collapses to a single line; where the dynamics settle onto a limit cycle it
# opens into a band, so a Hopf bifurcation appears as the point where the band
# starts to widen.
#
# It is the fishing-effort special case of `scanModel()`, which scans any aspect
# of a model and measures any quantity on the attractor it settles on.

#' Draw a bifurcation diagram over fishing effort
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Sweeps the fishing effort over a range of values and, for each value, runs
#' the full dynamics to their attractor and records the long-term range of a
#' summary quantity (biomass by default). The result is a bifurcation diagram
#' with fishing effort on the x-axis.
#'
#' This is [scanModel()] with [scanEffort()] as its setter, drawn with
#' `style = "envelope"`. Use `scanModel()` directly to scan something other than
#' the fishing effort, or to measure something other than biomass, yield or
#' spawning stock biomass.
#'
#' For each effort value the attractor is found in two stages. First
#' [projectToSteady()] runs the dynamics until they settle, stopping early once
#' it detects a stable steady state or a limit cycle and reporting which via its
#' `"convergence"` attribute. Then the settled state is projected forward once
#' more with [project()] over a sampling window — exactly one period for a limit
#' cycle, `t_sample` years for a run that settled on neither — and the minimum
#' and maximum of the chosen quantity over that window are taken as the
#' attractor envelope. A stable steady state does not move, so it is not sampled
#' at all and its minimum and maximum are equal.
#'
#' For a stable steady state the minimum and maximum coincide and the species is
#' drawn as a single line; once the dynamics settle onto a limit cycle the two
#' separate and the species is drawn as a shaded band between them. The onset of
#' the band therefore marks the Hopf bifurcation (see
#' `vignette("dynamic_stability")`).
#'
#' By default the sweep uses *continuation*: the projection at each effort value
#' starts from the attractor reached at the previous value rather than from the
#' original state. This shortens the transient and keeps the sweep on a single
#' branch of the attractor. Because it follows one branch, the diagram can look
#' different for an increasing versus a decreasing effort sequence if the model
#' has coexisting attractors (hysteresis); pass a decreasing `effort` vector to
#' trace the other direction.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param effort A numeric vector of fishing effort values for the x-axis. The
#'   same effort is applied to every gear (as in `project(effort = value)`).
#'   Default `seq(0, 2, length.out = 21)`.
#' @param species The species to include. By default all target species. A
#'   vector of species names or indices, as for other mizer plotting functions.
#' @param value The quantity for the y-axis, one of `"biomass"` (default),
#'   `"yield"` or `"ssb"`, computed with [getBiomass()], [getYield()] or
#'   [getSSB()] respectively.
#' @param style How the envelope is drawn, see [plot.MizerScan()]. The default
#'   `"envelope"` draws a line along the top and the bottom of the band.
#' @param t_max The maximum number of years to run the settling stage
#'   ([projectToSteady()]) at each effort value.
#' @param t_sample The length in years of the window over which the attractor is
#'   sampled when the run settled on neither a fixed point nor a limit cycle.
#'   Default `10`.
#' @param tol Convergence tolerance for the settling stage, passed to
#'   [projectToSteady()]. Tighter (smaller) values settle the stable branch more
#'   fully, giving cleaner single lines below the bifurcation at the cost of
#'   longer runs. Default `0.01`.
#' @param amplitude_tol The minimum relative biomass amplitude for a run to be
#'   classified as a limit cycle, passed to [projectToSteady()]. Default `0.01`.
#' @param extinction_threshold The relative reproduction collapse below which a
#'   species is treated as extinct, passed to [projectToSteady()]. Default
#'   `1e-6`.
#' @param continuation If `TRUE` (default) each settling run warm-starts from the
#'   attractor of the previous effort value.
#' @param return_data If `TRUE` the [MizerScan] object underlying the plot is
#'   returned instead of the plot. Default `FALSE`.
#' @param progress_bar If `TRUE` a text progress bar is shown while the effort
#'   values are swept. Defaults to `interactive()`.
#' @param log_y,log Whether to use a logarithmic y-axis, see [parsePlotLog()].
#' @param t_sample_default `r lifecycle::badge("deprecated")` Use `t_sample`
#'   instead.
#' @param t_save `r lifecycle::badge("deprecated")` Ignored: the attractor is
#'   now sampled at every time step.
#' @param ytrans `r lifecycle::badge("deprecated")` Use `log_y` instead.
#' @param ... Further arguments are passed on to [scanModel()].
#' @return A ggplot2 object, or, if `return_data = TRUE`, the [MizerScan] object
#'   holding the data. That object inherits from `data.frame` and has the
#'   columns `Effort`, the measured quantity, `Species`, `ymin`, `ymax`, `type`,
#'   `settled`, `period` and `residual`.
#' @seealso [scanModel()], [getStability()], [projectToSteady()], [plotBiomass()]
#' @family plotting functions
#' @family scan functions
#' @export
#' @examples
#' \donttest{
#' plotBifurcation(NS_params, effort = seq(0, 2, length.out = 11))
#' }
plotBifurcation <- function(params,
                            effort = seq(0, 2, length.out = 21),
                            species = NULL,
                            value = c("biomass", "yield", "ssb"),
                            style = "envelope",
                            t_max = 100,
                            t_sample = 10,
                            tol = 0.01,
                            amplitude_tol = 0.01,
                            extinction_threshold = 1e-6,
                            continuation = TRUE,
                            return_data = FALSE,
                            progress_bar = interactive(),
                            log_y = TRUE,
                            log = NULL,
                            t_sample_default = lifecycle::deprecated(),
                            t_save = lifecycle::deprecated(),
                            ytrans = lifecycle::deprecated(),
                            ...) {
    if (lifecycle::is_present(t_sample_default)) {
        lifecycle::deprecate_soft("3.2.1.9005",
                                  "plotBifurcation(t_sample_default)",
                                  "plotBifurcation(t_sample)")
        if (missing(t_sample)) t_sample <- t_sample_default
    }
    if (lifecycle::is_present(t_save)) {
        lifecycle::deprecate_soft(
            "3.2.1.9005", "plotBifurcation(t_save)",
            details = paste("The attractor is now sampled at every time step,",
                            "so `t_save` is ignored."))
    }
    if (lifecycle::is_present(ytrans)) {
        lifecycle::deprecate_soft("3.2.1.9005", "plotBifurcation(ytrans)",
                                  "plotBifurcation(log_y)")
        log_y <- identical(ytrans, "log10")
    }
    value <- match.arg(value)
    assert_that(is.numeric(effort), length(effort) >= 2)
    # Resolved here rather than left to `scanModel()`, which selects against
    # the series a value function returned and so knows nothing about
    # background species or about a logical or numeric species argument.
    params <- validParams(params)
    species <- valid_species_arg(params, species, error_on_empty = TRUE)

    scan <- scanModel(params, scan_values = effort,
                      set_func = scanEffort(),
                      value_func = switch(value,
                                          biomass = getBiomass,
                                          yield   = getYield,
                                          ssb     = getSSB),
                      species = species,
                      scan_name = "Effort",
                      t_max = t_max, t_sample = t_sample, tol = tol,
                      amplitude_tol = amplitude_tol,
                      extinction_threshold = extinction_threshold,
                      continuation = continuation,
                      progress_bar = progress_bar,
                      ...)

    if (return_data) return(scan)

    plot(scan, style = style, log_y = log_y, log = log)
}
