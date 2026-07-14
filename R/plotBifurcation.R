# Bifurcation diagram over fishing effort -------------------------------------
#
# `plotBifurcation()` sweeps the fishing effort, follows the attractor of the
# full dynamics at each effort value, and plots the long-term range (envelope)
# of a chosen summary quantity. Where the steady state is stable the envelope
# collapses to a single line; where the dynamics settle onto a limit cycle it
# opens into a band, so a Hopf bifurcation appears as the point where the band
# starts to widen.

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
#' For each effort value the attractor is found in two stages. First
#' [projectToSteady()] runs the dynamics until they settle, stopping early once
#' it detects a stable steady state or a limit cycle and reporting which via its
#' `"convergence"` attribute. Then the settled state is projected forward once
#' more with [project()] over a short sampling window — one full period for a
#' limit cycle, or a few years otherwise — and the minimum and maximum of the
#' chosen quantity over that window are taken as the attractor envelope.
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
#' @param t_max The maximum number of years to run the settling stage
#'   ([projectToSteady()]) at each effort value.
#' @param t_sample The length in years of the window over which the settled
#'   attractor is sampled to measure the envelope. If `NULL` (default) it is
#'   chosen automatically: one full period for a detected limit cycle, or
#'   `t_sample_default` years otherwise.
#' @param t_sample_default The sampling window used when no limit-cycle period is
#'   available (a stable or non-converged run). Default `10`.
#' @param t_save The interval at which the sampling window is saved, controlling
#'   how finely the cycle envelope is resolved. Default `0.25`.
#' @param continuation If `TRUE` (default) each settling run warm-starts from the
#'   attractor of the previous effort value.
#' @param return_data If `TRUE` the data frame underlying the plot is returned
#'   instead of the plot. Default `FALSE`.
#' @param progress_bar If `TRUE` (default) a text progress bar is shown while the
#'   effort values are swept.
#' @param ytrans Transformation for the y-axis, `"log10"` (default) or
#'   `"identity"`.
#' @return A ggplot2 object, or a data frame with columns `Effort`, `Species`,
#'   `ymin`, `ymax` and `type` (the attractor type reported by
#'   [projectToSteady()]) if `return_data = TRUE`.
#' @seealso [getStability()], [projectToSteady()], [plotBiomass()]
#' @family plotting functions
#' @export
#' @examples
#' \donttest{
#' plotBifurcation(NS_params, effort = seq(0, 2, length.out = 11))
#' }
plotBifurcation <- function(params,
                            effort = seq(0, 2, length.out = 21),
                            species = NULL,
                            value = c("biomass", "yield", "ssb"),
                            t_max = 100,
                            t_sample = NULL,
                            t_sample_default = 10,
                            t_save = 0.25,
                            continuation = TRUE,
                            return_data = FALSE,
                            progress_bar = TRUE,
                            ytrans = "log10") {
    params <- validParams(params)
    value <- match.arg(value)
    assert_that(is.numeric(effort), length(effort) >= 2,
                is.number(t_max), t_max > 0,
                is.null(t_sample) || is.number(t_sample),
                is.number(t_sample_default), t_sample_default > 0)
    species <- valid_species_arg(params, species, error_on_empty = TRUE)

    value_func <- switch(value,
                         biomass = getBiomass,
                         yield   = getYield,
                         ssb     = getSSB)
    ylab <- switch(value,
                   biomass = "Biomass [g]",
                   yield   = "Yield [g/year]",
                   ssb     = "Spawning stock biomass [g]")

    if (isTRUE(progress_bar)) {
        pb <- utils::txtProgressBar(min = 0, max = length(effort), style = 3)
        on.exit(close(pb), add = TRUE)
    }

    p_run <- params
    rows <- vector("list", length(effort))
    for (i in seq_along(effort)) {
        # Stage 1: settle onto the attractor, stopping early once a steady state
        # or limit cycle is detected.
        settled <- projectToSteady(p_run, effort = effort[i], t_max = t_max,
                                   return_sim = FALSE, progress_bar = FALSE,
                                   info_level = 0)
        conv <- attr(settled, "convergence")

        # Stage 2: sample the settled attractor to measure the envelope. Starting
        # from the settled state means the window contains no transient.
        window <- if (!is.null(t_sample)) {
            t_sample
        } else if (identical(conv$type, "cycle")) {
            ceiling(conv$period)
        } else {
            t_sample_default
        }
        sim <- project(settled, effort = effort[i], t_max = window,
                       t_save = t_save, progress_bar = FALSE)
        vals <- value_func(sim)[, species, drop = FALSE]

        rows[[i]] <- data.frame(
            Effort  = effort[i],
            Species = species,
            ymin    = apply(vals, 2, min),
            ymax    = apply(vals, 2, max),
            type    = conv$type,
            row.names = NULL,
            stringsAsFactors = FALSE
        )

        # Warm-start the next effort value from the attractor just reached.
        if (continuation) p_run <- settled
        if (isTRUE(progress_bar)) utils::setTxtProgressBar(pb, i)
    }
    plot_dat <- do.call(rbind, rows)

    if (return_data) return(plot_dat)

    species_levels <- intersect(names(params@linecolour),
                                unique(plot_dat$Species))
    plot_dat$Species <- factor(plot_dat$Species, levels = species_levels)
    linecolour <- params@linecolour[species_levels]

    ybreaks <- if (ytrans == "log10") log_breaks() else waiver()

    p <- ggplot(plot_dat, aes(x = .data$Effort, group = .data$Species)) +
        geom_ribbon(aes(ymin = .data$ymin, ymax = .data$ymax,
                        fill = .data$Species),
                    alpha = 0.25, colour = NA) +
        geom_line(aes(y = .data$ymax, colour = .data$Species)) +
        geom_line(aes(y = .data$ymin, colour = .data$Species)) +
        scale_colour_manual(values = linecolour) +
        scale_fill_manual(values = linecolour) +
        scale_y_continuous(trans = ytrans, breaks = ybreaks,
                           labels = prettyNum, name = ylab) +
        scale_x_continuous(name = "Fishing effort")

    make_mizer_plot(p, c("Species", "Effort", "ymax", "ymin"))
}
