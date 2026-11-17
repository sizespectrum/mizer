# Size spectra, their comparisons and relative differences.

#' Convert a weight-based spectrum to a length-based spectrum
#'
#' A density with respect to weight is converted to a density with respect to
#' length with the Jacobian `dw/dl = b * w / l`. A density with respect to
#' logarithmic weight is instead converted with
#' `d log(w) / d log(l) = b`. This is the interface used by the `power`-based
#' spectrum plots; arrays carry their density measure explicitly and use
#' [convert_plot_density_axis()] instead.
#'
#' @inheritParams convert_plot_density_axis
#' @param power The power of weight multiplying the number density.
#' @param per_log_size Whether the spectrum is a density with respect to
#'   logarithmic size rather than with respect to size. Defaults to
#'   `power == 2`, the only power for which this used to be the case.
#' @return The plotting data with both its size coordinate and spectrum values
#'   expressed for the requested axis.
#' @keywords internal
convert_plot_spectrum_axis <- function(plot_dat, params, size_axis, power,
                                       per_log_size = power == 2,
                                       species_col = "Species",
                                       value_col = 2) {
    convert_plot_density_axis(plot_dat, params, size_axis,
                              density_wrt = spectrum_density_wrt(per_log_size),
                              species_col = species_col,
                              value_col = value_col)
}

#' The density measure of a `power`-based spectrum
#'
#' @param per_log_size Whether the spectrum is a density with respect to
#'   logarithmic size.
#' @return `"log_w"` or `"w"`.
#' @keywords internal
spectrum_density_wrt <- function(per_log_size) {
    if (isTRUE(per_log_size)) "log_w" else "w"
}

#' Plot abundance and biomass spectra
#'
#' `plotSpectra()` plots either a number density or a biomass density, either
#' with respect to size or with respect to logarithmic size. Those two choices
#' are made with the `biomass` and `per_log_size` arguments. When called with a
#' [MizerSim] object, the abundance is averaged over the specified time range
#' (a single value for the time range can be used to plot a single time step).
#' When called with a [MizerParams] object the initial abundance is plotted.
#' With `size_axis = "l"`, densities are converted from per unit weight to per
#' unit length; densities with respect to logarithmic size are instead
#' converted between logarithmic weight and logarithmic length intervals.
#'
#' The plotted quantity is the number density multiplied by `w^power`, where
#' the power is the sum of the two choices above: a biomass density carries one
#' factor of the weight and a density with respect to logarithmic size carries
#' another:
#'
#' | | `per_log_size = FALSE` | `per_log_size = TRUE` |
#' |---|---|---|
#' | `biomass = FALSE` | `power = 0` | `power = 1` |
#' | `biomass = TRUE`  | `power = 1` | `power = 2` |
#'
#' The `power` argument can still be given instead, and is the only way to ask
#' for a power that is not the sum of the two flags. But note that `power` on
#' its own does not distinguish the two entries with `power = 1`: it is taken
#' to mean the biomass density with respect to weight, which is what determines
#' the y-axis label and the Jacobian used for a length axis. Supplying `power`
#' together with a flag that contradicts it is an error.
#'
#' The `log_x` argument only controls how the size axis is displayed; it does
#' not change the density on the y-axis. In particular, showing weight on a
#' logarithmic axis does not by itself convert a density per unit weight into a
#' density per logarithmic weight interval. That choice is made with
#' `per_log_size`, and the conversion from weight to length then uses the
#' logarithmic Jacobian, irrespective of the value of `log_x`.
#'
#' `plotlySpectra()` is the interactive plotly version. To compare spectra from
#' two objects use [plotSpectra2()]. To show relative differences use
#' [plotSpectraRelative()].
#'
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the w axis. Use NA for the default: the lower default is
#'   `min(params@w) / 100` when `resource = TRUE` (to show some resource below
#'   the fish grid) or `min(params@w)` when `resource = FALSE`; the upper
#'   default is `max(params@w_full)`. These limits only set the visible window;
#'   they do not filter the returned data or change any calculation.
#' @param llim A numeric vector of length two providing lower and upper limits
#'   for the length axis when `size_axis = "l"`. Use `NA` to auto-scale to the
#'   data range. These limits only set the visible window.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use `NA` to auto-scale to the data range. These limits
#'   only set the visible window, without filtering the returned data.
#' @param biomass Whether to plot the biomass density (`TRUE`, the default) or
#'   the number density (`FALSE`).
#' @param per_log_size Whether to plot the density with respect to logarithmic
#'   size (`TRUE`) or with respect to size (`FALSE`, the default).
#' @param power The abundance is plotted as the number density times the weight
#'   raised to `power`. An alternative to the `biomass` and `per_log_size`
#'   arguments, with which it must agree if they are given as well; see
#'   Details. The default is `power = 1`, the biomass density.
#' @param total A boolean value that determines whether the total is plotted as
#'   well. The total is the total of everything the object holds — every
#'   species and the resource — whatever is drawn, so it does not move when
#'   `species`, `resource` or `background` change. Default is FALSE.
#' @param resource A boolean value that determines whether resource is included.
#'   Default is TRUE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted by
#'   being plotted with thicker lines.
#' @param log_x If `TRUE` (default), use a log10 x-axis.
#' @param log_y If `TRUE` (default), use a log10 y-axis.
#' @param log Character string specifying which axes should use log10 scales,
#'   in the same form as the base [plot()] argument. For example, `"x"`,
#'   `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.
#' @param size_axis Whether to plot size as weight (`"w"`, default) or length
#'   (`"l"`), using the allometric weight-length relationship. Spectrum
#'   densities and their units are transformed to match the chosen axis.
#' @param return_data A boolean value that determines whether the formatted data
#' used for the plot is returned instead of the plot itself. Default value is FALSE
#' @param ... Further arguments used by only some of the methods:
#'
#'   **For `MizerSim` methods:**
#'   * `time_range`: The time range (either a vector of values, a vector
#'     of min and max time, or a single value) to average the abundances
#'     over. Default is the final time step.
#'   * `geometric_mean`: `r lifecycle::badge("experimental")` If `TRUE`
#'     then the average of the abundances over the time range is a geometric
#'     mean instead of the default arithmetic mean.
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'w' (or 'l' if `size_axis = "l"`), 'value',
#'   'Species', 'Legend' is returned. `plotlySpectra()` returns a plotly object.
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @examples
#' \donttest{
#' params <-  NS_params
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotSpectra(sim)
#' plotSpectra(sim, wlim = c(1e-6, NA))
#' plotSpectra(sim, time_range = 10:20)
#' plotSpectra(sim, time_range = 10:20, biomass = FALSE)
#' plotSpectra(sim, species = c("Cod", "Herring"), per_log_size = TRUE)
#' plotSpectra(sim, species = c("Cod", "Herring"), size_axis = "l")
#'
#' # Returning the data frame
#' fr <- plotSpectra(sim, return_data = TRUE)
#' str(fr)
#' }
#' @rdname plotSpectra
#' @export
plotSpectra <- function(object, species = NULL,
                        wlim = c(NA, NA), llim = c(NA, NA), ylim = c(NA, NA),
                        power = NULL, biomass = NULL, per_log_size = NULL,
                        total = FALSE,
                        resource = TRUE, background = TRUE, highlight = NULL,
                        log_x = TRUE, log_y = TRUE, log = NULL,
                        size_axis = c("w", "l"), return_data = FALSE, ...) {
    UseMethod("plotSpectra")
}

#' @rdname plotSpectra
#' @usage NULL
#' @export
plotSpectra.MizerSim <- function(object, species = NULL,
                        wlim = c(NA, NA), llim = c(NA, NA),
                        ylim = c(NA, NA),
                        power = NULL, biomass = NULL, per_log_size = NULL,
                        total = FALSE, resource = TRUE,
                        background = TRUE,
                        highlight = NULL, log_x = TRUE, log_y = TRUE,
                        log = NULL, size_axis = c("w", "l"),
                        return_data = FALSE,
                        time_range,
                        geometric_mean = FALSE, ...) {
    spectrum <- resolve_spectrum_power(power, biomass, per_log_size)
    size_axis <- plot_size_axis(size_axis)
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    assert_that(is.flag(total), is.flag(resource),
                is.flag(background),
                length(wlim) == 2,
                length(llim) == 2,
                length(ylim) == 2)
    species <- valid_species_arg(object, species)
    if (length(species) == 0 && !total && !resource) {
        stop("There is nothing to plot as no valid species have been selected.")
    }
    abundance <- sim_mean_abundance(object, time_range, geometric_mean)
    plot_spectra(object@params, n = abundance$n, n_pp = abundance$n_pp,
                 species = species, wlim = wlim, ylim = ylim,
                 llim = llim,
                 power = spectrum$power, biomass = spectrum$biomass,
                 per_log_size = spectrum$per_log_size,
                 total = total, resource = resource,
                 background = background, highlight = highlight,
                 log_x = log_x, log_y = log_y,
                 size_axis = size_axis,
                 return_data = return_data)
}

#' @rdname plotSpectra
#' @usage NULL
#' @export
plotSpectra.MizerParams <- function(object, species = NULL,
                        wlim = c(NA, NA), llim = c(NA, NA),
                        ylim = c(NA, NA),
                        power = NULL, biomass = NULL, per_log_size = NULL,
                        total = FALSE, resource = TRUE,
                        background = TRUE,
                        highlight = NULL, log_x = TRUE, log_y = TRUE,
                        log = NULL, size_axis = c("w", "l"),
                        return_data = FALSE, ...) {
    spectrum <- resolve_spectrum_power(power, biomass, per_log_size)
    size_axis <- plot_size_axis(size_axis)
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    assert_that(is.flag(total), is.flag(resource),
                is.flag(background),
                length(wlim) == 2,
                length(llim) == 2,
                length(ylim) == 2)
    species <- valid_species_arg(object, species)
    if (length(species) == 0 && !total && !resource) {
        stop("There is nothing to plot as no valid species have been selected.")
    }
    plot_spectra(object, n = object@initial_n,
                 n_pp = object@initial_n_pp,
                 species = species, wlim = wlim, ylim = ylim,
                 llim = llim,
                 power = spectrum$power, biomass = spectrum$biomass,
                 per_log_size = spectrum$per_log_size,
                 total = total, resource = resource,
                 background = background, highlight = highlight,
                 log_x = log_x, log_y = log_y,
                 size_axis = size_axis,
                 return_data = return_data)
}

#' Build the size-spectrum plot
#'
#' Internal worker shared by the [plotSpectra()] methods. It assembles the
#' plotting data frame from the species and resource abundances, optionally
#' converts to a length axis, and either
#' returns the data or draws the plot via [plotDataFrame()].
#'
#' @param params A MizerParams object.
#' @param n Array of species abundances (species by size).
#' @param n_pp Vector of resource abundance.
#' @param species The species to be plotted.
#' @param wlim,llim Numeric vectors of length two giving the weight and length
#'   limits.
#' @param ylim Numeric vector of length two giving the y-axis limits.
#' @param power The abundance is multiplied by weight raised to this power.
#' @param biomass Whether the resulting quantity is a biomass density rather
#'   than a number density. Used for the y-axis label.
#' @param per_log_size Whether the resulting quantity is a density with respect
#'   to logarithmic size. Used for the y-axis label and for the Jacobian of the
#'   conversion to a length axis.
#' @param total Whether to include the total community abundance.
#' @param resource Whether to include the resource spectrum.
#' @param background Whether to include background species.
#' @param highlight Name or vector of names of species to be highlighted.
#' @param log_x,log_y Logical flags for log10 axes.
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @param return_data If `TRUE`, return the plotting data frame instead of the
#'   plot.
#' @return A `mizer_plot` (ggplot2) object, or the plotting data frame if
#'   `return_data = TRUE`.
#' @keywords internal
plot_spectra <- function(params, n, n_pp,
                         species, wlim, llim, ylim, power,
                         biomass = power >= 1, per_log_size = power == 2,
                         total, resource, background,
                         highlight, log_x, log_y, size_axis, return_data) {
    params <- validParams(params)
    size_axis <- plot_size_axis(size_axis)
    # The spectrum is a density: each plotted value N_j w^power is a point on the
    # continuous N(w) w^power curve, which (for the cell-average N_j) it touches
    # at the geometric bin centre. Under second-order bin-averaging we therefore
    # evaluate *both* the w^power weight and the plotted location at the bin
    # centre w* = w sqrt(beta) (issue #383). The default keeps the grid nodes, so
    # existing spectrum plots are unchanged.
    second_order <- isTRUE(params@second_order_w[["bin_average"]])
    w_grid      <- if (second_order) bin_midpoints(params) else params@w
    w_full_grid <- if (second_order) bin_midpoints(params, w_full = TRUE) else
        params@w_full
    # Default size limits follow the plotted grid (the bin centres under
    # second order) so the centre shift does not clip the top/bottom bins. With
    # the default (first-order) nodes these reduce to the previous limits.
    if (is.na(wlim[1])) {
        wlim[1] <- if (resource) min(w_grid) / 100 else min(w_grid)
    }
    if (is.na(wlim[2])) {
        wlim[2] <- max(w_full_grid)
    }

    y_label <- spectra_y_label(power, size_axis, biomass = biomass,
                               per_log_size = per_log_size)
    series <- abundance_series(params, n, n_pp, power, w_grid, w_full_grid,
                               background = background)
    plot_dat <- series[series$Species %in%
                           shown_series(params, species, resource, background),
                       , drop = FALSE]
    # Species first, then the resource, then the background.
    block <- ifelse(plot_dat$Legend == "Background", 2,
                    ifelse(plot_dat$Species == "Resource", 1, 0))
    plot_dat <- plot_dat[order(block), , drop = FALSE]
    plot_dat <- plot_dat[plot_dat$value > 0, ]
    plot_dat <- convert_plot_spectrum_axis(plot_dat, params, size_axis, power,
                                           per_log_size = per_log_size)
    if (total) {
        # The total is the total of everything the model holds: every species,
        # whether or not it was selected for display, and the resource, whether
        # or not it is drawn. It is summed over the series as they are plotted,
        # which on a length axis means at equal length rather than at equal
        # weight. Zeros are kept here, unlike in the plotted data: they are
        # part of each series and dropping them would let the interpolation
        # bridge a gap that is really empty.
        total_dat <- series
        total_dat$Legend <- "Total"
        total_dat <- convert_plot_spectrum_axis(total_dat, params, size_axis,
                                                power,
                                                per_log_size = per_log_size)
        total_dat <- add_total_line(total_dat,
                                    x_var = plot_size_x_var(size_axis),
                                    value_col = "value")
        total_dat <- total_dat[total_dat$Species == "Total" &
                                   total_dat$value > 0, ]
        plot_dat <- rbind(plot_dat, total_dat[, names(plot_dat), drop = FALSE])
    }
    # Retain the existing negligible-density floor, independently of limits.
    plot_dat <- plot_dat[plot_dat$value > 1e-20, ]
    names(plot_dat)[2] <- y_label
    if (return_data) return(plot_dat)

    plotDataFrame(plot_dat, params, xlab = plot_size_xlab(size_axis),
                  ylab = y_label,
                  xtrans = if (log_x) "log10" else "identity",
                  ytrans = if (log_y) "log10" else "identity",
                  xlim = plot_size_xlim(wlim, size_axis, llim), ylim = ylim,
                  highlight = highlight, legend_var = "Legend")
}

# The abundances of a simulation averaged over `time_range` (by default the
# final saved time), arithmetically or geometrically, as plotSpectra() and
# plotCDF() show them.
sim_mean_abundance <- function(sim, time_range, geometric_mean = FALSE) {
    if (missing(time_range)) {
        time_range <- max(as.numeric(dimnames(sim@n)$time))
    }
    time_elements <- get_time_elements(sim, time_range)
    mean_fn <- if (geometric_mean) function(x) exp(mean(log(x))) else mean
    list(n = apply(sim@n[time_elements, , , drop = FALSE], c(2, 3), mean_fn),
         n_pp = apply(sim@n_pp[time_elements, , drop = FALSE], 2, mean_fn))
}

# One row for every size bin of every species and of the resource, holding the
# abundance times `weight^power`. The rows are placed at `x` (`x_full` for the
# resource) and the power of weight is evaluated at `weight` (`weight_full`):
# a spectrum uses one grid for both, a cumulative distribution places rows on
# the left bin edges but weights them where the scheme does. Background
# species carry the "Background" legend when `background = TRUE`.
abundance_series <- function(params, n, n_pp, power, x, x_full,
                             weight = x, weight_full = x_full,
                             background = TRUE) {
    species <- as.character(params@species_params$species)
    series <- rbind(
        data.frame(w = rep(x, each = length(species)),
                   value = c(sweep(unclass(n), 2, weight^power, "*")),
                   Species = rep(species, times = length(x))),
        data.frame(w = x_full, value = c(n_pp) * weight_full^power,
                   Species = "Resource")
    )
    series$Legend <- series$Species
    if (background) {
        background_species <- species[params@species_params$is_background]
        series$Legend[series$Species %in% background_species] <- "Background"
    }
    series
}

# The series a spectrum or cumulative plot draws: the selected species, every
# background species when `background = TRUE` (once, under the "Background"
# legend, even when it was also selected by name) and the resource when
# `resource = TRUE`.
shown_series <- function(params, species, resource, background) {
    shown <- valid_species_arg(params, species)
    if (background) {
        species_params <- params@species_params
        shown <- union(shown, as.character(
            species_params$species[species_params$is_background]))
    }
    if (resource) shown <- c(shown, "Resource")
    shown
}

#' Compare abundance and biomass spectra from two objects
#'
#' `plotSpectra2()` compares the abundance spectra from two `MizerParams` or
#' `MizerSim` objects in a single plot. Colours identify species or groups and
#' linetype identifies the object.
#'
#' `plotlySpectra2()` is the interactive plotly version.
#'
#' @param object1 First `MizerParams` or `MizerSim` object.
#' @param object2 Second `MizerParams` or `MizerSim` object.
#' @param name1,name2 Labels for the two objects, used in the linetype legend.
#' @inheritParams plotSpectra
#' @param ... Additional arguments passed to [plotSpectra()] for preparing the
#'   spectra data, for example `time_range` or `geometric_mean` for `MizerSim`
#'   objects.
#'
#' @return A ggplot2 object. `plotlySpectra2()` returns a plotly object.
#'
#' @family plotting functions
#' @seealso [plotting_functions], [plotSpectra()], [plotSpectraRelative()]
#' @examples
#' \donttest{
#' sim1 <- project(NS_params, t_max = 10, progress_bar = FALSE)
#' sim2 <- project(NS_params, effort = 0.5, t_max = 10, progress_bar = FALSE)
#' plotSpectra2(sim1, sim2, "Original", "Effort = 0.5")
#' }
#' @export
plotSpectra2 <- function(object1, object2, name1 = "First", name2 = "Second",
                         species = NULL,
                         wlim = c(NA, NA), llim = c(NA, NA),
                         ylim = c(NA, NA),
                         power = NULL, biomass = NULL, per_log_size = NULL,
                         total = FALSE, resource = TRUE,
                         background = TRUE,
                         highlight = NULL,
                         log_x = TRUE, log_y = TRUE,
                         log = NULL, size_axis = c("w", "l"), ...) {
    size_axis <- plot_size_axis(size_axis)
    spectrum <- resolve_spectrum_power(power, biomass, per_log_size)
    power <- spectrum$power
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y
    assert_that(length(wlim) == 2, length(llim) == 2, length(ylim) == 2)

    # The size axis is converted by `plotSpectra()` rather than afterwards, so
    # that the total it forms is the total on the axis being plotted. A total
    # summed at equal weight cannot be converted to a length axis afterwards:
    # it is a sum over species that each convert differently. It also means
    # `ylim` can be handed straight over, because the values it is applied to
    # are then the values it describes.
    sf1 <- plotSpectra(object1, species = species, wlim = wlim, ylim = ylim,
                       power = power, total = total, resource = resource,
                       background = background, size_axis = size_axis,
                       return_data = TRUE, ...)
    sf2 <- plotSpectra(object2, species = species, wlim = wlim, ylim = ylim,
                       power = power, total = total, resource = resource,
                       background = background, size_axis = size_axis,
                       return_data = TRUE, ...)
    params <- if (is(object1, "MizerSim")) object1@params else object1

    plotComparisonDataFrame(sf1, sf2, params,
                            name1 = name1, name2 = name2,
                            xlab = plot_size_xlab(size_axis),
                            ylab = spectra_y_label(
                                power, size_axis,
                                biomass = spectrum$biomass,
                                per_log_size = spectrum$per_log_size),
                            xtrans = if (log_x) "log10" else "identity",
                            ytrans = if (log_y) "log10" else "identity",
                            xlim = plot_size_xlim(wlim, size_axis, llim),
                            ylim = ylim, highlight = highlight,
                            legend_var = "Legend")
}

#' Resolve the power of weight multiplying a spectrum
#'
#' The quantity plotted by [plotSpectra()] is the number density multiplied by
#' `w^power`. That power is the sum of two independent choices: whether the
#' quantity is a biomass density (a factor of `w`) or a number density, and
#' whether it is a density with respect to logarithmic size (another factor of
#' `w`) or with respect to size. The two choices are what determine the y-axis
#' label and the Jacobian used when converting to a length axis, and they are
#' not recoverable from `power` alone: `power = 1` is both the biomass density
#' with respect to weight and the number density with respect to logarithmic
#' weight.
#'
#' Users therefore express the choice with the `biomass` and `per_log_size`
#' flags. The `power` argument remains available, both for backwards
#' compatibility and as an escape hatch for powers that are not the sum of two
#' flags. Each argument is `NULL` when it was not supplied by the user.
#'
#' @param power The power of weight multiplying the number density, or `NULL`.
#' @param biomass Whether to plot a biomass density rather than a number
#'   density, or `NULL`. The default is `TRUE`.
#' @param per_log_size Whether to plot a density with respect to logarithmic
#'   size rather than with respect to size, or `NULL`. The default is `FALSE`.
#' @return A list with entries `power`, `biomass` and `per_log_size`.
#' @keywords internal
resolve_spectrum_power <- function(power = NULL, biomass = NULL,
                                   per_log_size = NULL) {
    if (!is.null(biomass)) assert_that(is.flag(biomass), noNA(biomass))
    if (!is.null(per_log_size)) {
        assert_that(is.flag(per_log_size), noNA(per_log_size))
    }
    if (is.null(power)) {
        biomass <- biomass %||% TRUE
        per_log_size <- per_log_size %||% FALSE
        return(list(power = as.numeric(biomass + per_log_size),
                    biomass = biomass, per_log_size = per_log_size))
    }
    assert_that(is.number(power), noNA(power))
    if (is.null(biomass) && is.null(per_log_size)) {
        # Only `power` was given, so infer the flags from it. This reproduces
        # the interpretation that mizer has always placed on `power`.
        return(list(power = power,
                    biomass = power >= 1, per_log_size = power == 2))
    }
    # `power` was given together with at least one flag. Deduce the other flag
    # from it and then insist that the three agree, rather than silently
    # letting one of them win (issue #501).
    if (is.null(biomass)) biomass <- (power - per_log_size) == 1
    if (is.null(per_log_size)) per_log_size <- (power - biomass) == 1
    if (!isTRUE(all.equal(power, biomass + per_log_size))) {
        stop("You gave `power = ", power, "` but `biomass = ", biomass,
             "` and `per_log_size = ", per_log_size, "` imply `power = ",
             biomass + per_log_size, "`. The power of weight multiplying the ",
             "number density is `biomass + per_log_size`, so supply either ",
             "`power` or the flags, but not contradictory values of both.")
    }
    list(power = power, biomass = biomass, per_log_size = per_log_size)
}

#' Y-axis label for a size-spectrum plot
#'
#' @param power The power of weight that the abundance was multiplied by.
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @param biomass Whether the quantity is a biomass density rather than a
#'   number density. Defaults to `power >= 1`.
#' @param per_log_size Whether the quantity is a density with respect to
#'   logarithmic size. Defaults to `power == 2`.
#' @return A character string for the y-axis label.
#' @keywords internal
spectra_y_label <- function(power, size_axis = "w",
                            biomass = power >= 1,
                            per_log_size = power == 2) {
    size_axis <- plot_size_axis(size_axis)
    if (power %in% c(0, 1, 2)) {
        length_axis <- identical(size_axis, "l")
        if (biomass) {
            if (per_log_size) return("Biomass density [g]")
            return(if (length_axis) "Biomass density [g/cm]" else
                "Biomass density")
        }
        if (per_log_size) {
            return(if (length_axis) "Number density in log length" else
                "Number density in log weight")
        }
        return(if (length_axis) "Number density [1/cm]" else
            "Number density [1/g]")
    }
    paste0("Number density * w^", power)
}

#' @rdname plotSpectra2
#' @usage NULL
#' @export
plotlySpectra2 <- function(object1, object2, name1 = "First",
                           name2 = "Second",
                           species = NULL,
                           wlim = c(NA, NA), llim = c(NA, NA),
                           ylim = c(NA, NA),
                           power = NULL, biomass = NULL, per_log_size = NULL,
                           total = FALSE, resource = TRUE,
                           background = TRUE,
                           highlight = NULL,
                           log_x = TRUE, log_y = TRUE, log = NULL,
                           size_axis = c("w", "l"), ...) {
    plotly_forward("plotSpectra2")
}

#' Plot relative difference between abundance spectra
#'
#' `plotSpectraRelative()` plots the difference between the spectra relative to
#' their average. If we denote the number density from the first object as
#' \eqn{N_1(w)} and that from the second object as \eqn{N_2(w)}, then this
#' plot shows
#' \deqn{2 (N_2(w) - N_1(w)) / (N_2(w) + N_1(w)).}
#' Note that it does not matter whether the relative difference is calculated
#' for number density, biomass density, or biomass density in log weight,
#' because the factors of \eqn{w} by which the densities differ cancel out in
#' the relative difference.
#'
#' `plotlySpectraRelative()` is the interactive plotly version.
#'
#' @param object1 First `MizerParams` or `MizerSim` object.
#' @param object2 Second `MizerParams` or `MizerSim` object.
#' @inheritParams plotSpectra
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis (the relative difference). Use `NA` to refer to the
#'   existing minimum or maximum.
#' @param ... Additional arguments passed to [plotSpectra()] for preparing the
#'   spectra data, for example `time_range` or `geometric_mean` for `MizerSim`
#'   objects.
#'
#' @return A ggplot2 object. `plotlySpectraRelative()` returns a plotly object.
#'
#' @family plotting functions
#' @seealso [plotting_functions], [plotSpectra()], [plotSpectra2()]
#' @examples
#' \donttest{
#' sim1 <- project(NS_params, t_max = 10, progress_bar = FALSE)
#' sim2 <- project(NS_params, effort = 0.5, t_max = 10, progress_bar = FALSE)
#' plotSpectraRelative(sim1, sim2)
#' }
#' @export
plotSpectraRelative <- function(object1, object2,
                                species = NULL,
                                wlim = c(NA, NA), llim = c(NA, NA),
                                ylim = c(NA, NA),
                                total = FALSE, resource = TRUE,
                                background = TRUE,
                                highlight = NULL,
                                log_x = TRUE,
                                size_axis = c("w", "l"), ...) {
    size_axis <- plot_size_axis(size_axis)
    assert_that(length(wlim) == 2, length(llim) == 2, length(ylim) == 2)

    # As in `plotSpectra2()`, the conversion is left to `plotSpectra()` so that
    # the total is summed on the axis being plotted.
    sf1 <- plotSpectra(object1, species = species, wlim = wlim,
                       power = 1, total = total, resource = resource,
                       background = background, size_axis = size_axis,
                       return_data = TRUE, ...)
    sf2 <- plotSpectra(object2, species = species, wlim = wlim,
                       power = 1, total = total, resource = resource,
                       background = background, size_axis = size_axis,
                       return_data = TRUE, ...)
    params <- if (is(object1, "MizerSim")) object1@params else object1
    params <- validParams(params)

    # A size axis, so the two spectra are interpolated onto a common grid.
    # On a weight axis they already share one and nothing is approximated; on a
    # length axis they need not, because the two models can convert weight to
    # length differently.
    plotRelativeDataFrame(sf1, sf2, params,
                          xlab = plot_size_xlab(size_axis),
                          xtrans = if (log_x) "log10" else "identity",
                          xlim = plot_size_xlim(wlim, size_axis, llim),
                          ylim = ylim,
                          highlight = highlight,
                          legend_var = "Legend", interpolate = TRUE)
}

#' @rdname plotSpectraRelative
#' @usage NULL
#' @export
plotlySpectraRelative <- function(object1, object2,
                                  species = NULL,
                                  wlim = c(NA, NA), llim = c(NA, NA),
                                  ylim = c(NA, NA),
                                  total = FALSE, resource = TRUE,
                                  background = TRUE,
                                  highlight = NULL,
                                  log_x = TRUE,
                                  size_axis = c("w", "l"), ...) {
    plotly_forward("plotSpectraRelative",
                   tooltip = plot_size_tooltip(size_axis, before = "Legend",
                                               after = "rel_diff"))
}

#' @rdname plotSpectra
#' @usage NULL
#' @export
plotlySpectra <- function(object, species = NULL,
                        time_range, geometric_mean = FALSE,
                        wlim = c(NA, NA), llim = c(NA, NA),
                        ylim = c(NA, NA),
                        power = NULL, biomass = NULL, per_log_size = NULL,
                        total = FALSE, resource = TRUE,
                        background = TRUE,
                        highlight = NULL, log_x = TRUE, log_y = TRUE,
                        log = NULL,
                        size_axis = c("w", "l"), ...) {
    plotly_forward("plotSpectra")
}
