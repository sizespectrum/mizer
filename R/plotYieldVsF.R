# Yield against the fishing mortality on one species -----------------------
#
# `plotYieldVsF()` is the one scan common enough in fisheries work to be worth
# its own name: the yield of a single stock plotted against the fishing
# mortality exerted on it, with the fishing on every other species left alone.
# The fishing mortality at the peak of that curve is F_MSY.
#
# It is `scanModel()` with `scanFishingMortality()` as its setter and
# `getYield()` as the quantity measured. Everything that makes the picture
# readable — the axis labels, the `F_MSY` reference line, the sweep outwards
# from the current fishing mortality — is already carried by the setter or by
# `scanModel()`, so nothing here is reimplemented.

#' Plot the yield of a species against the fishing mortality on it
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Varies the fishing mortality on one species over a range of values, leaving
#' the fishing on every other species unchanged, and plots the long-term yield
#' of that species against it. The fishing mortality at which the yield is
#' largest is \eqn{F_{MSY}}, and is marked on the plot by default.
#'
#' This is [scanModel()] with [scanFishingMortality()] as its setter and
#' [getYield()] as the quantity it measures. Use `scanModel()` directly to vary
#' something other than the fishing mortality on a single species, to measure
#' something other than the yield, or to follow more than one species at once.
#'
#' At each fishing mortality the model is projected until it settles, and what
#' is plotted depends on what it settled on. At a fixed point the yield is read
#' straight off the settled state. On a limit cycle it is averaged over exactly
#' one period, and the band around the line shows the range the yield covers
#' over that cycle, so an oscillation is displayed rather than silently averaged
#' away. Fishing mortalities at which the model settled on neither are marked
#' with a cross and should not be relied on; raise `t_max` for those.
#'
#' The scan starts from the fishing mortality the model currently sits at and
#' works outwards in both directions, each arm warm-starting from the attractor
#' reached at the previous value.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param species The name of the species whose fishing mortality is varied.
#'   Only one species at a time.
#' @param F_range A numeric vector of fishing mortalities for the x-axis. If
#'   missing it is built as `seq(F_min, F_max, length.out = no_steps)`.
#' @param F_min,F_max,no_steps Used to build `F_range` when that is missing.
#' @param gear The name of the gear whose fishing mortality on the species is
#'   varied. Only needed when several gears catch the species; if NULL
#'   (default), the fishing mortality from all of them is replaced. See
#'   [scanFishingMortality()].
#' @param style How the range covered on a limit cycle is drawn, see
#'   [plot.MizerScan()]. The default `"ribbon"` draws the average as a line
#'   inside the band.
#' @param mark_max Whether to mark the fishing mortality at which the yield is
#'   largest, which is \eqn{F_{MSY}}. Default TRUE.
#' @param reference_lines Whether to draw the `F_MSY` species parameter, if the
#'   species has one, as a vertical line. See [plot.MizerScan()].
#' @param log_y,log Whether to use a logarithmic y-axis, see [parsePlotLog()].
#'   Unlike most mizer plots this defaults to FALSE, because the yield is
#'   exactly zero at zero fishing mortality and a logarithmic axis would drop
#'   the point that anchors the curve.
#' @param return_data If TRUE the [MizerScan] object underlying the plot is
#'   returned instead of the plot. Default FALSE.
#' @param progress_bar If TRUE a text progress bar is shown while the fishing
#'   mortalities are swept. Defaults to `interactive()`.
#' @param ... Further arguments are passed on to [scanModel()].
#'
#' @return A ggplot2 object, or, if `return_data = TRUE`, the [MizerScan] object
#'   holding the data. The fishing mortality giving the largest yield is
#'   available from that object as `attr(scan, "at_max")`.
#' @seealso [scanModel()], [scanFishingMortality()], [plot.MizerScan()],
#'   [getYield()]
#' @family plotting functions
#' @family scan functions
#' @concept scan
#' @export
#' @examples
#' \donttest{
#' plotYieldVsF(NS_params, "Cod", F_max = 1.5, no_steps = 8)
#'
#' # The fishing mortality that maximises the yield
#' scan <- plotYieldVsF(NS_params, "Cod", F_max = 1.5, no_steps = 8,
#'                      return_data = TRUE)
#' attr(scan, "at_max")
#' }
plotYieldVsF <- function(params, species, F_range,
                         F_min = 0, F_max = 1.5, no_steps = 16,
                         gear = NULL,
                         style = "ribbon",
                         mark_max = TRUE,
                         reference_lines = TRUE,
                         log_y = FALSE,
                         log = NULL,
                         return_data = FALSE,
                         progress_bar = interactive(),
                         ...) {
    params <- validParams(params)
    # Resolved here rather than left to the setter, so that a mistyped species
    # or a species argument selecting more than one is caught before any
    # projection has run.
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) != 1) {
        stop("You can only plot the yield against the fishing mortality for ",
             "one species at a time, but ", length(species),
             " species were selected. Use `scanModel()` to follow several ",
             "species through the same scan.")
    }
    if (missing(F_range)) {
        assert_that(is.number(F_min), is.number(F_max), F_max > F_min,
                    is.number(no_steps), no_steps >= 2)
        F_range <- seq(F_min, F_max, length.out = no_steps)
    }
    assert_that(is.numeric(F_range), length(F_range) >= 2)

    scan <- scanModel(params, scan_values = F_range,
                      set_func = scanFishingMortality(species, gear = gear),
                      value_func = getYield,
                      species = species,
                      # Sweep outwards from where the model already sits, so
                      # that neither arm has to travel through the other.
                      current_scan_value = "auto",
                      progress_bar = progress_bar,
                      ...)

    if (return_data) return(scan)

    plot(scan, style = style, mark_max = mark_max,
         reference_lines = reference_lines, log_y = log_y, log = log)
}
