# Biomass time series and the interactive compatibility wrapper.

#' Plot the biomass of species through time
#'
#' After running a projection, the biomass of each species can be plotted
#' against time. The biomass is calculated within user defined size limits
#' (see [getBiomass()]).
#'
#' @param object An object of class \linkS4class{MizerSim}
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param tlim A numeric vector of length two providing lower and upper limits
#'   for the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no limit at
#'   that end. Default is `c(NA, NA)`.
#' @param y_ticks The approximate number of ticks desired on the y axis.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use `NA` to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param total A boolean value that determines whether the total biomass from
#'   all species is plotted as well. Default is FALSE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param log_x If `TRUE`, use a log10 x-axis. Default is `FALSE`.
#' @param log_y If `TRUE`, use a log10 y-axis. Default is `TRUE`.
#' @param log Character string specifying which axes should use log10 scales,
#'   in the same form as the base [plot()] argument. For example, `"x"`,
#'   `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.
#'   For backward compatibility, `TRUE` and `FALSE` are interpreted as setting
#'   only `log_y`.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default is FALSE.
#' @param use_cutoff If TRUE, the `biomass_cutoff` column in the species
#'   parameters is used as the minimum weight for each species.
#' @param ... Arguments setting the size range over which the biomass is
#'   calculated (see [getBiomass()]):
#'   \describe{
#'     \item{`min_w`}{Smallest weight in size range. Defaults to smallest weight
#'       in the model.}
#'     \item{`max_w`}{Largest weight in size range. Defaults to largest weight
#'       in the model.}
#'     \item{`min_l`}{Smallest length in size range. If supplied, this takes
#'       precedence over `min_w`.}
#'     \item{`max_l`}{Largest length in size range. If supplied, this takes
#'       precedence over `max_w`.}
#'   }
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'Year', 'Biomass', 'Species', 'Legend' is
#'   returned.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions], [getBiomass()]
#' @examples
#' \donttest{
#' plotBiomass(NS_sim)
#' plotBiomass(NS_sim, species = c("Sandeel", "Herring"), total = TRUE)
#' plotBiomass(NS_sim, tlim = c(1980, 1990))
#'
#' # Returning the data frame
#' fr <- plotBiomass(NS_sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotBiomass
#' @export
plotBiomass <- function(object, species = NULL, tlim = c(NA, NA),
                        y_ticks = 6, ylim = c(NA, NA), total = FALSE,
                        background = TRUE, highlight = NULL, log = NULL,
                        return_data = FALSE, log_x = FALSE, log_y = TRUE,
                        use_cutoff = FALSE, ...) {
    UseMethod("plotBiomass")
}

#' @rdname plotBiomass
#' @usage NULL
#' @export
plotBiomass.MizerSim <- function(object, species = NULL,
                        tlim = c(NA, NA),
                        y_ticks = 6, ylim = c(NA, NA),
                        total = FALSE, background = TRUE,
                        highlight = NULL,
                        log = NULL, return_data = FALSE,
                        log_x = FALSE, log_y = TRUE,
                        use_cutoff = FALSE,
                        start_time = lifecycle::deprecated(),
                        end_time = lifecycle::deprecated(),
                        min_w = min(object@params@w),
                        max_w = max(object@params@w),
                        min_l = NULL, max_l = NULL, ...) {
    if (lifecycle::is_present(start_time)) {
        lifecycle::deprecate_warn("2.6.0", "plotBiomass(start_time)", "plotBiomass(tlim)")
        tlim[[1]] <- start_time
    }
    if (lifecycle::is_present(end_time)) {
        lifecycle::deprecate_warn("2.6.0", "plotBiomass(end_time)", "plotBiomass(tlim)")
        tlim[[2]] <- end_time
    }
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    bm <- getBiomass(object, use_cutoff = use_cutoff,
                     min_w = min_w, max_w = max_w,
                     min_l = min_l, max_l = max_l)
    plot(bm, species = species,
         tlim = tlim,
         y_ticks = y_ticks, ylim = ylim,
         total = total, background = background,
         highlight = highlight,
         log_x = log_axes$log_x, log_y = log_axes$log_y,
         return_data = return_data)
}

#' @rdname plotBiomass
#' @usage NULL
#' @export
plotlyBiomass <- function(object,
             species = NULL,
             tlim = c(NA, NA),
             start_time = lifecycle::deprecated(),
             end_time = lifecycle::deprecated(),
             y_ticks = 6,
             ylim = c(NA, NA),
             total = FALSE,
             background = TRUE,
             highlight = NULL,
             log = NULL, log_x = FALSE, log_y = TRUE,
             use_cutoff = FALSE,
             min_w = min(object@params@w),
             max_w = max(object@params@w),
             min_l = NULL,
             max_l = NULL) {
    if (lifecycle::is_present(start_time)) {
        lifecycle::deprecate_warn("2.6.0", "plotlyBiomass(start_time)", "plotlyBiomass(tlim)")
        tlim[[1]] <- start_time
    }
    if (lifecycle::is_present(end_time)) {
        lifecycle::deprecate_warn("2.6.0", "plotlyBiomass(end_time)", "plotlyBiomass(tlim)")
        tlim[[2]] <- end_time
    }
    plotly_forward("plotBiomass", tooltip = c("Species", "Year", "Biomass"),
                   exclude = c("start_time", "end_time"))
}
