# Plotting functions ----

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2019 Gustav Delius and Richard Southwell.
# Distributed under the GPL 3 or later
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

#' Description of the plotting functions
#'
#' Mizer provides a range of plotting functions for visualising the results
#' of running a simulation, stored in a MizerSim object, or the initial state
#' stored in a MizerParams object.
#'
#' The quickest way to make a standard plot is often to call `plot()` directly.
#' mizer provides `plot()` methods for `MizerSim` and `MizerParams` objects, and
#' also for the array classes returned by many summary and rate functions:
#'
#' * `plot(<MizerSim>)` produces a five-panel summary plot with
#'   [plotFeedingLevel()], [plotBiomass()], [plotPredMort()], [plotFMort()] and
#'   [plotSpectra()].
#'
#' * `plot(<MizerParams>)` produces a three-panel summary plot with
#'   [plotFeedingLevel()], [plotPredMort()] and [plotSpectra()] from the initial
#'   state.
#'
#' * `plot(<ArrayTimeBySpecies>)` plots any time-by-species array, such as
#'   those returned by [getBiomass()], [getSSB()], [getYield()] and [getN()] on
#'   a `MizerSim`, as lines of value against time.
#'
#' * `plot(<ArraySpeciesBySize>)` plots any species-by-size array, such as
#'   those returned by [getEncounter()], [getFeedingLevel()], [getPredMort()],
#'   [getFMort()] and [search_vol()], as lines of value against body size.
#'
#' * `plot(<ArrayTimeBySpeciesBySize>)` plots a time slice from a
#'   time-by-species-by-size array, such as those returned by [getFMort()] or
#'   [getPredMort()] on a `MizerSim`.
#'
#' * `plot(<ArrayResourceBySize>)` plots a resource array, such as
#'   those returned by [getResourceMort()], against body size.
#'
#' * `plot(<ArrayTimeByResourceBySize>)` plots a time slice from a
#'   time-by-resource array, such as that returned by [NResource()] on a
#'   `MizerSim`, against body size.
#'
#' The same array objects can be passed to [plotHover()] to produce
#' hover-enabled plotly versions, for example `plotHover(getBiomass(sim))` or
#' `plotHover(getEncounter(params))`. To add another compatible array to an
#' existing ggplot, use [addPlot()]. To compare two compatible mizer arrays
#' directly, use [plot2()]. To plot cumulative distributions over body size,
#' use [plotCDF()]. To visualise how spectra or rates change through time, use
#' [animate()] on a `MizerSim` or an
#' `ArrayTimeBySpeciesBySize` object.
#'
#' The named plotting functions give more specialised control. This table shows
#' the available named plotting functions.
#' \tabular{ll}{
#'   Plot \tab Description \cr
#'   [plotBiomass()] \tab Plots the total biomass of each species through time. A time range to be plotted can be specified. The size range of the community can be specified in the same way as for [getBiomass()]. \cr
#'   [plotYield()] \tab Plots the total yield of each species across all fishing gears against time. \cr
#'   [plotYieldGear()] \tab Plots the total yield of each species by gear against time. \cr
#'   [plotSpectra()] \tab Plots the abundance (biomass or numbers) spectra of each species and the background community. It is possible to specify a minimum size which is useful for truncating the plot. \cr
#'   [plotCDF()] \tab Plots cumulative distributions of abundance or biomass over size. \cr
#'   [plotCDF2()] \tab Compares cumulative distributions from two simulations or parameter objects in one plot. \cr
#'   [plotSpectra2()] \tab Compares the spectra from two simulations or parameter objects in one plot. \cr
#'   [plotFeedingLevel()] \tab Plots the feeding level of each species against size. \cr
#'   [plotPredMort()] \tab Plots the predation mortality of each species against size. \cr
#'   [plotFMort()] \tab Plots the total fishing mortality of each species against size. \cr
#'   [plotGrowthCurves()] \tab Plots the size as a function of age. \cr
#'   [plotDiet()] \tab Plots the diet composition at size for a given predator species. \cr
#'   [plotBiomassObservedVsModel()] \tab Compares observed biomass with model biomass. \cr
#'   [plotYieldObservedVsModel()] \tab Compares observed yield with model yield. \cr
#'   [animate()] \tab Animates spectra or rate arrays through time. The older [animateSpectra()] name is retained as an alias. \cr
#' }
#'
#' The static plotting functions use ggplot2 and return a ggplot object. This
#' means that you can manipulate the plot further after its creation using the
#' ggplot grammar of graphics. The named high-level plot functions have plotly
#' counterparts, for example [plotlyBiomass()] or [plotlySpectra()], for
#' interactive exploration. Generic and compositional plotting APIs, such as
#' [plot()], [plot2()], [plotRelative()] and [addPlot()], do not have separate
#' plotly wrappers. Use [plotHover()] on the ggplot object they return.
#'
#' While most plot functions take their data from a MizerSim object, some of
#' those that make plots representing data at a single time can also take their
#' data from the initial values in a MizerParams object.
#'
#' Where plots show results for species, the line colour and line type for each
#' species are specified by the `linecolour` and `linetype` slots in
#' the MizerParams object. These were either taken from a default palette
#' hard-coded into [emptyParams()] or they were specified by the user
#' in the species parameters dataframe used to set up the MizerParams object.
#' The `linecolour` and `linetype` slots hold named vectors, named by
#' the species. They can be overwritten by the user at any time.
#'
#' Most plots allow the user to select to show only a subset of species,
#' specified as a vector in the `species` argument to the plot function.
#'
#' The ordering of the species in the legend is the same as the ordering in
#' the species parameter data frame.
#'
#' @seealso [summary_functions], [indicator_functions]
#' @family plotting functions
#' @name plotting_functions
#' @examples
#' \donttest{
#' sim <- NS_sim
#'
#' # Generic plot methods
#' plot(sim)
#' plot(getBiomass(sim), species = c("Cod", "Herring"))
#' plotHover(getBiomass(sim))
#'
#' # Named plot functions
#' plotFeedingLevel(sim)
#'
#' # Plotting only a subset of species
#' plotFeedingLevel(sim, species = c("Cod", "Herring"))
#'
#' # Adding another compatible array to an existing plot
#' p <- plot(getBiomass(sim), species = "Cod")
#' addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")
#'
#' # Specifying new colours and linetypes for some species
#' sim@params@linetype["Cod"] <- "dashed"
#' sim@params@linecolour["Cod"] <- "red"
#' plotFeedingLevel(sim, species = c("Cod", "Herring"))
#'
#' # Manipulating the plot
#' library(ggplot2)
#' p <- plotFeedingLevel(sim)
#' p <- p + geom_hline(aes(yintercept = 0.7))
#' p <- p + theme_bw()
#' p
#' }
NULL

#### plot ####

# Arguments of a summary plot that only the spectrum panel understands.
# `plot()` on a MizerParams or MizerSim object passes its `...` on to every
# panel it draws, but these describe the plotted spectrum and mean nothing to a
# rate. They must not reach the other panels: `per_log_size` is deliberately an
# error on an array that is not a density (see `check_per_log_size()`), so
# broadcasting it made `plot(sim, per_log_size = TRUE)` fail on the feeding
# level panel while the equivalent `power = 2` went through.
spectra_only_args <- c("power", "biomass", "per_log_size", "resource")

shared_plot_args <- function(dots) {
    if (length(dots) == 0) return(dots)
    nms <- names(dots) %||% rep("", length(dots))
    dots[!(nms %in% spectra_only_args)]
}

#' Summary plot for `MizerSim` objects
#'
#' After running a projection, produces 5 plots in the same window: feeding
#' level, abundance spectra, predation mortality and fishing mortality of each
#' species by size; and biomass of each species through time. This method just
#' puts the plots generated by [plotBiomass()], [plotFeedingLevel()],
#' [plotSpectra()], [plotPredMort()] and [plotFMort()] all in one window.
#'
#' @param x An object of class \linkS4class{MizerSim}
#' @param ...  Arguments passed on to the individual plotting functions
#'   [plotBiomass()], [plotFeedingLevel()], [plotSpectra()], [plotPredMort()]
#'   and [plotFMort()]. The arguments that describe the plotted spectrum
#'   (`power`, `biomass`, `per_log_size` and `resource`) go to [plotSpectra()]
#'   only, because they mean nothing to a rate.
#' @return A viewport object
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @rdname plotMizerSim
#' @name plotMizerSim
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plot(sim)
#' }
plot.MizerSim <- function(x, ...) {
    shared <- shared_plot_args(list(...))
    p1 <- do.call(plotFeedingLevel, c(list(x), shared))
    p2 <- plotSpectra(x, ...)
    p3 <- do.call(plotBiomass, c(list(x, y_ticks = 3), shared))
    p4 <- do.call(plotPredMort, c(list(x), shared))
    p5 <- do.call(plotFMort, c(list(x), shared))
    grid::grid.newpage()
    glayout <- grid::grid.layout(3, 2) # widths and heights arguments
    vp <- grid::viewport(layout = glayout)
    grid::pushViewport(vp)
    vplayout <- function(x, y) {
        grid::viewport(layout.pos.row = x, layout.pos.col = y)
    }
    print(p1 + theme(legend.position = "none"), vp = vplayout(1, 1))
    print(p3 + theme(legend.position = "none"), vp = vplayout(1, 2))
    print(p4 + theme(legend.position = "none"), vp = vplayout(2, 1))
    print(p5 + theme(legend.position = "none"), vp = vplayout(2, 2))
    print(p2 + theme(legend.position = "right",
                     legend.key.size = unit(0.1, "cm")),
          vp = vplayout(3, 1:2))
}

#' Summary plot for `MizerParams` objects
#'
#' Produces 3 plots in the same window: abundance spectra, feeding
#' level and predation mortality of each species against size. This method just
#' puts the plots generated by [plotFeedingLevel()], [plotPredMort()] and
#' [plotSpectra()] all in one window.
#'
#' @param x An object of class \linkS4class{MizerParams}
#' @param ...  Arguments passed on to the individual plotting functions
#'   [plotFeedingLevel()], [plotSpectra()] and [plotPredMort()]. The arguments
#'   that describe the plotted spectrum (`power`, `biomass`, `per_log_size` and
#'   `resource`) go to [plotSpectra()] only, because they mean nothing to a
#'   rate.
#' @return A viewport object
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @name plotMizerParams
#' @examples
#' \donttest{
#' params <-  NS_params
#' plot(params)
#' }
plot.MizerParams <- function(x, ...) {
    params <- validParams(x)
    shared <- shared_plot_args(list(...))
    p11 <- do.call(plotFeedingLevel, c(list(params), shared))
    p2 <- plotSpectra(params, ...)
    p12 <- do.call(plotPredMort, c(list(params), shared))
    grid::grid.newpage()
    glayout <- grid::grid.layout(2, 2) # widths and heights arguments
    vp <- grid::viewport(layout = glayout)
    grid::pushViewport(vp)
    vplayout <- function(x, y) {
        grid::viewport(layout.pos.row = x, layout.pos.col = y)
    }
    print(p11 + theme(legend.position = "none"), vp = vplayout(1, 1))
    print(p12 + theme(legend.position = "none"), vp = vplayout(1, 2))
    print(p2 + theme(legend.position = "right",
                     legend.key.size = unit(0.1, "cm")),
          vp = vplayout(2, 1:2))
}
