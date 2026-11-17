# Cumulative distributions and explicit integration bounds.

#' Plot cumulative abundance or biomass distributions
#'
#' `plotCDF()` plots the cumulative distribution over body size from small to
#' large sizes. The complete number density is multiplied by `w^power` and
#' integrated over the range chosen by `min_w`, `max_w`, `min_l` and `max_l`.
#' Plot limits only change the visible window; they do not change the
#' cumulative values, normalisation, totals or `return_data`.
#' With `normalise = TRUE`, each curve is divided by its final value so that it
#' ends at 1.
#'
#' Unlike for [plotSpectra()], the only choice that matters here is `biomass`:
#' whether to accumulate numbers or biomass. Whether a density is expressed
#' with respect to size or with respect to logarithmic size makes no difference
#' to its integral, because the change of variable cancels the factor of the
#' weight, and so `plotCDF()` does not accept `per_log_size = TRUE`.
#'
#' `plotlyCDF()` is the interactive plotly version. To compare cumulative
#' distributions from two objects, use [plotCDF2()].
#'
#' @param object An object of class \linkS4class{MizerSim} or
#'   \linkS4class{MizerParams}.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param min_w,max_w Single non-negative numbers giving the weight range in
#'   grams to integrate over. Defaults are 0 and `Inf`, covering the complete
#'   grid.
#'   A boundary inside a bin includes only the overlapping fraction of its
#'   width. Density is weighted at the left bin boundary in first-order models
#'   and at the geometric midpoint when bin averaging is enabled.
#'
#'   The bounds are exact, which is not yet true of the `min_w` and `max_w` of
#'   [getBiomass()] and the other summary functions: in a first-order model
#'   they include a whole bin when its left edge lies in the range, so
#'   `getBiomass(params, max_w = w)` with `w` on a bin edge still includes the
#'   bin that starts at `w`, while a cumulative curve with `max_w = w` ends at
#'   `w`. With bin averaging the two also weight each bin slightly differently.
#' @param min_l,max_l Optional single non-negative numbers giving the length
#'   range in centimetres to integrate over. Each overrides the corresponding
#'   weight bound, using each species' own length-weight relationship, including
#'   the resource's equivalent-sphere convention. Independent of `size_axis`.
#' @param wlim A numeric vector of length two giving the visible weight-axis
#'   limits. Use `NA` to auto-scale. Does not select the integration range.
#' @param llim A numeric vector of length two providing lower and upper limits
#'   for the length axis when `size_axis = "l"`. Use `NA` to auto-scale to the
#'   data range. These limits only set the visible window.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use `NA` to auto-scale to the data range. These limits
#'   only set the visible window, without filtering the returned data.
#' @param biomass Whether to plot the cumulative biomass (`TRUE`, the default)
#'   or the cumulative abundance (`FALSE`).
#' @param per_log_size Only `FALSE` (the default) is accepted; see Details.
#' @param power The number density is multiplied by the weight raised to
#'   `power` before being integrated. An alternative to the `biomass`
#'   argument, with which it must agree if that is given as well. The default
#'   is `power = 1`, the cumulative biomass.
#' @param total A boolean value that determines whether the total is plotted as
#'   well. The total is the total of everything the object holds — every
#'   species and the resource — whatever is drawn, so it does not move when
#'   `species`, `resource` or `background` change. Default is FALSE.
#' @param resource A boolean value that determines whether resource is included.
#'   Default is FALSE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted by
#'   being plotted with thicker lines.
#' @param normalise If `TRUE` (default), plot the cumulative proportion. If
#'   `FALSE`, plot the cumulative abundance, biomass, or other unnormalised
#'   integral. A curve with zero mass in the integration range has `NA`
#'   normalised values (and zero unnormalised values).
#' @param log_x If `TRUE` (default), use a log10 x-axis.
#' @param log_y If `TRUE`, use a log10 y-axis. Default is `FALSE`.
#' @param log Character string specifying which axes should use a log10 scale,
#'   in the same form as the base [plot()] argument. If supplied, this overrides
#'   `log_x` and `log_y`.
#' @param size_axis Whether to plot size as weight (`"w"`, default) or length
#'   (`"l"`), using the allometric weight-length relationship.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default is FALSE.
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
#'   'Species', 'Legend' is returned. `plotlyCDF()` returns a plotly object.
#' @export
#' @family plotting functions
#' @seealso [plotSpectra()], [plotCDF2()]
#' @examples
#' \donttest{
#' plotCDF(NS_params, species = c("Cod", "Herring"))
#' plotCDF(NS_sim, power = 0, normalise = FALSE)
#' }
plotCDF <- function(object, species = NULL,
                    wlim = c(NA, NA), llim = c(NA, NA), ylim = c(NA, NA),
                    power = NULL, biomass = NULL, per_log_size = NULL,
                    total = FALSE,
                    resource = FALSE, background = TRUE, highlight = NULL,
                    normalise = TRUE, log_x = TRUE, log_y = FALSE, log = NULL,
                    size_axis = c("w", "l"), return_data = FALSE,
                    min_w = 0, max_w = Inf, min_l = NULL, max_l = NULL, ...) {
    UseMethod("plotCDF")
}

#' @rdname plotCDF
#' @usage NULL
#' @export
plotCDF.MizerSim <- function(object, species = NULL,
                             wlim = c(NA, NA), llim = c(NA, NA),
                             ylim = c(NA, NA),
                             power = NULL, biomass = NULL,
                             per_log_size = NULL,
                             total = FALSE, resource = FALSE,
                             background = TRUE,
                             highlight = NULL, normalise = TRUE,
                             log_x = TRUE, log_y = FALSE, log = NULL,
                             size_axis = c("w", "l"),
                             return_data = FALSE,
                             time_range,
                             geometric_mean = FALSE,
                             min_w = 0, max_w = Inf,
                             min_l = NULL, max_l = NULL, ...) {
    abundance <- sim_mean_abundance(object, time_range, geometric_mean)
    params <- object@params
    params@initial_n[] <- abundance$n
    params@initial_n_pp[] <- abundance$n_pp
    plotCDF(params, species = species, wlim = wlim, llim = llim, ylim = ylim,
            power = power, biomass = biomass, per_log_size = per_log_size,
            total = total, resource = resource, background = background,
            highlight = highlight, normalise = normalise,
            log_x = log_x, log_y = log_y, log = log, size_axis = size_axis,
            return_data = return_data, min_w = min_w, max_w = max_w,
            min_l = min_l, max_l = max_l, ...)
}

#' @rdname plotCDF
#' @usage NULL
#' @export
plotCDF.MizerParams <- function(object, species = NULL,
                                wlim = c(NA, NA), llim = c(NA, NA),
                                ylim = c(NA, NA),
                                power = NULL, biomass = NULL,
                                per_log_size = NULL,
                                total = FALSE, resource = FALSE,
                                background = TRUE,
                                highlight = NULL, normalise = TRUE,
                                log_x = TRUE, log_y = FALSE, log = NULL,
                                size_axis = c("w", "l"),
                                return_data = FALSE,
                                min_w = 0, max_w = Inf,
                                min_l = NULL, max_l = NULL, ...) {
    size_axis <- plot_size_axis(size_axis)
    power <- resolve_cdf_power(power, biomass, per_log_size)$power
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    assert_that(is.flag(total), is.flag(resource),
                is.flag(background), is.flag(normalise),
                length(wlim) == 2,
                length(llim) == 2,
                length(ylim) == 2)

    plot_cdf(object, power = power, normalise = normalise,
             log_x = log_axes$log_x, log_y = log_axes$log_y,
             wlim = wlim, ylim = ylim, llim = llim,
             highlight = highlight, size_axis = size_axis,
             return_data = return_data, species = species, total = total,
             resource = resource, background = background,
             min_w = min_w, max_w = max_w, min_l = min_l, max_l = max_l)
}

# Build CDFs from complete, unfiltered abundance data. A display selection must
# never decide which species or resource bins contribute to the community total.
plot_cdf <- function(params, power, normalise, log_x, log_y, wlim, llim,
                     ylim, highlight, size_axis, return_data, species,
                     total, resource, background,
                     min_w = 0, max_w = Inf, min_l = NULL, max_l = NULL) {
    params <- validParams(params)
    shown <- shown_series(params, species, resource, background)
    if (total) shown <- c(shown, "Total")
    # The density is weighted once, at the left bin edge or at the geometric
    # bin centre according to bin_average; the rows stay on the left edges so
    # that prepare_spectra_cdf_data() can find each bin.
    second_order <- isTRUE(params@second_order_w[["bin_average"]])
    series <- abundance_series(
        params, params@initial_n, params@initial_n_pp, power,
        x = params@w, x_full = params@w_full,
        weight = if (second_order) bin_midpoints(params) else params@w,
        weight_full = if (second_order) {
            bin_midpoints(params, w_full = TRUE)
        } else {
            params@w_full
        },
        background = background)
    # A total includes every contributor, whichever curves are drawn.
    if (!total) series <- series[series$Species %in% shown, , drop = FALSE]
    cdf_dat <- prepare_spectra_cdf_data(series, params, normalise = FALSE,
                                        min_w = min_w, max_w = max_w,
                                        min_l = min_l, max_l = max_l)
    cdf_dat <- convert_plot_size_axis(cdf_dat, params, size_axis)
    if (total) {
        cdf_dat <- add_total_line(cdf_dat, x_var = plot_size_x_var(size_axis),
                                  value_col = "value", cumulative = TRUE)
    }
    cdf_dat <- cdf_dat[cdf_dat$Species %in% shown, , drop = FALSE]
    if (nrow(cdf_dat) == 0) {
        stop("The CDF integration bounds do not overlap ",
             "the selected size grids.")
    }
    if (normalise) cdf_dat <- normalise_cdf(cdf_dat)
    names(cdf_dat)[names(cdf_dat) == "value"] <- cdf_y_label(power, normalise)
    if (return_data) return(cdf_dat)

    plotDataFrame(cdf_dat, params,
                  xlab = plot_size_xlab(size_axis),
                  ylab = cdf_y_label(power, normalise),
                  xtrans = if (log_x) "log10" else "identity",
                  ytrans = if (log_y) "log10" else "identity",
                  xlim = plot_size_xlim(wlim, size_axis, llim), ylim = ylim,
                  highlight = highlight, legend_var = "Legend")
}

# Convert scalar integration bounds to weights for each contributing species.
# Length bounds override the corresponding weight bounds, just as they do in
# get_size_range_array(). The resource uses its own length-weight relationship.
cdf_weight_bounds <- function(params, species, min_w, max_w, min_l, max_l) {
    bounds <- list(min_w = min_w, max_w = max_w, min_l = min_l, max_l = max_l)
    for (name in names(bounds)) {
        value <- bounds[[name]]
        if (is.null(value) && name %in% c("min_l", "max_l")) next
        if (!is.numeric(value) || length(value) != 1 || is.na(value) ||
                value < 0 || (startsWith(name, "min") && !is.finite(value))) {
            stop("`", name, "` must be a single non-negative number",
                 if (startsWith(name, "min")) " that is finite." else ".")
        }
    }
    lower <- rep(min_w, length(species))
    upper <- rep(max_w, length(species))
    if (!is.null(min_l) || !is.null(max_l)) {
        lw <- plot_length_params(species, params)
        if (anyNA(lw) || any(!is.finite(as.matrix(lw))) ||
                any(lw$a <= 0 | lw$b <= 0)) {
            stop("CDF length bounds need positive finite ",
                 "length-weight parameters.")
        }
        if (!is.null(min_l)) lower <- lw$a * min_l^lw$b
        if (!is.null(max_l)) upper <- lw$a * max_l^lw$b
    }
    if (any(lower >= upper)) {
        stop("Each CDF lower integration bound must be less than ",
             "its upper bound.")
    }
    list(lower = lower, upper = upper)
}

# The density has already been weighted once, at the left boundary or geometric
# midpoint according to bin_average. Integrate that bin density over the exact
# overlap with the requested interval: never average or apply the mask again.
# Cumulative values belong on the upper edge of that overlap, including a
# partial last bin. Full bins retain the existing CDF quadrature under both
# schemes, and zero-density bins retain their places in the cumulative curve.
prepare_spectra_cdf_data <- function(plot_dat, params, normalise = TRUE,
                                     min_w = 0, max_w = Inf,
                                     min_l = NULL, max_l = NULL) {
    plot_dat <- plot_dat[order(plot_dat$Species, plot_dat$w), ]
    bounds <- cdf_weight_bounds(params, plot_dat$Species,
                                min_w, max_w, min_l, max_l)
    widths <- spectra_bin_width(plot_dat$w, params)
    lower <- pmax(plot_dat$w, bounds$lower)
    upper <- pmin(plot_dat$w + widths, bounds$upper)
    overlap <- pmax(0, upper - lower)
    full_bin <- lower == plot_dat$w & upper == plot_dat$w + widths
    overlap[full_bin] <- widths[full_bin]
    plot_dat$value <- plot_dat$value * overlap
    plot_dat$w <- upper
    plot_dat <- plot_dat[overlap > 0, , drop = FALSE]
    if (nrow(plot_dat) == 0) return(plot_dat)
    plot_dat$value <- stats::ave(plot_dat$value, plot_dat$Species,
                                 FUN = cumsum)
    if (normalise) plot_dat <- normalise_cdf(plot_dat)
    plot_dat
}

normalise_cdf <- function(plot_dat) {
    totals <- stats::ave(plot_dat$value, plot_dat$Species, FUN = max)
    plot_dat$value <- ifelse(totals > 0, plot_dat$value / totals, NA_real_)
    plot_dat
}

#' Look up the size-bin widths for spectra data
#'
#' Matches the weights in the plotting data to the model's full size grid and
#' returns the corresponding bin widths.
#'
#' @param w Numeric vector of weights.
#' @param params A MizerParams object.
#' @return A numeric vector of bin widths, one for each entry of `w`.
#' @keywords internal
spectra_bin_width <- function(w, params) {
    idx <- match(w, params@w_full)
    if (anyNA(idx)) {
        missing <- which(is.na(idx))
        for (i in missing) {
            idx[i] <- which.min(abs(params@w_full - w[i]))
        }
        if (!isTRUE(all.equal(w, params@w_full[idx], scale = 1))) {
            stop("Could not determine size-bin widths for the spectra data.")
        }
    }
    params@dw_full[idx]
}

#' Y-axis label for a cumulative-distribution plot
#'
#' @param power The power of weight that the abundance was multiplied by.
#' @param normalise Whether the cumulative distribution is normalised.
#' @return A character string for the y-axis label.
#' @keywords internal
cdf_y_label <- function(power, normalise) {
    if (normalise) {
        if (power == 0) {
            return("Cumulative proportion of abundance")
        }
        if (power == 1) {
            return("Cumulative proportion of biomass")
        }
        return("Cumulative proportion")
    }
    if (power == 0) {
        return("Cumulative abundance")
    }
    if (power == 1) {
        return("Cumulative biomass [g]")
    }
    paste0("Cumulative number density * w^", power)
}

#' Compare cumulative abundance or biomass distributions from two objects
#'
#' `plotCDF2()` compares cumulative distributions from two `MizerParams` or
#' `MizerSim` objects in a single plot. Colours identify species or groups and
#' linetype identifies the object.
#'
#' `plotlyCDF2()` is the interactive plotly version.
#'
#' @param object1 First `MizerParams` or `MizerSim` object.
#' @param object2 Second `MizerParams` or `MizerSim` object.
#' @param name1,name2 Labels for the two objects, used in the linetype legend.
#' @inheritParams plotCDF
#' @param ... Additional arguments passed to [plotCDF()] for preparing the
#'   cumulative distribution data, for example `time_range` or `geometric_mean`
#'   for `MizerSim` objects.
#'
#' @return A ggplot2 object. `plotlyCDF2()` returns a plotly object.
#'
#' @family plotting functions
#' @seealso [plotSpectra()], [plotCDF()]
#' @examples
#' \donttest{
#' sim1 <- project(NS_params, t_max = 10, progress_bar = FALSE)
#' sim2 <- project(NS_params, effort = 0.5, t_max = 10, progress_bar = FALSE)
#' plotCDF2(sim1, sim2, "Original", "Effort = 0.5")
#' }
#' @export
plotCDF2 <- function(object1, object2, name1 = "First", name2 = "Second",
                     species = NULL,
                     wlim = c(NA, NA), llim = c(NA, NA),
                     ylim = c(NA, NA),
                     power = NULL, biomass = NULL, per_log_size = NULL,
                     total = FALSE, resource = FALSE,
                     background = TRUE,
                     highlight = NULL,
                     normalise = TRUE,
                     log_x = TRUE, log_y = FALSE,
                     log = NULL, size_axis = c("w", "l"),
                     min_w = 0, max_w = Inf,
                     min_l = NULL, max_l = NULL, ...) {
    size_axis <- plot_size_axis(size_axis)
    power <- resolve_cdf_power(power, biomass, per_log_size)$power
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    assert_that(is.flag(normalise),
                length(wlim) == 2, length(llim) == 2, length(ylim) == 2)

    cf1 <- plotCDF(object1, species = species,
                   power = power, total = total, resource = resource,
                   background = background, normalise = normalise,
                   size_axis = size_axis, return_data = TRUE,
                   min_w = min_w, max_w = max_w,
                   min_l = min_l, max_l = max_l, ...)
    cf2 <- plotCDF(object2, species = species,
                   power = power, total = total, resource = resource,
                   background = background, normalise = normalise,
                   size_axis = size_axis, return_data = TRUE,
                   min_w = min_w, max_w = max_w,
                   min_l = min_l, max_l = max_l, ...)
    params <- if (is(object1, "MizerSim")) object1@params else object1

    plotComparisonDataFrame(cf1, cf2, params,
                            name1 = name1, name2 = name2,
                            xlab = plot_size_xlab(size_axis),
                            ylab = cdf_y_label(power, normalise),
                            xtrans = if (log_axes$log_x) "log10" else "identity",
                            ytrans = if (log_axes$log_y) "log10" else "identity",
                            xlim = plot_size_xlim(wlim, size_axis, llim),
                            ylim = ylim, highlight = highlight,
                            legend_var = "Legend")
}

#' Resolve the power of weight for a cumulative distribution
#'
#' As [resolve_spectrum_power()], except that `per_log_size = TRUE` is
#' rejected: integrating a density with respect to logarithmic size gives the
#' same cumulative quantity as integrating the corresponding density with
#' respect to size, so the flag would be meaningless here.
#'
#' @inheritParams resolve_spectrum_power
#' @return A list with entries `power`, `biomass` and `per_log_size`.
#' @keywords internal
resolve_cdf_power <- function(power = NULL, biomass = NULL,
                              per_log_size = NULL) {
    if (isTRUE(per_log_size)) {
        stop("A cumulative distribution does not depend on whether the ",
             "density is taken with respect to size or with respect to ",
             "logarithmic size, so `per_log_size = TRUE` is not allowed here. ",
             "Use `biomass` to choose between the cumulative abundance and ",
             "the cumulative biomass.")
    }
    spectrum <- resolve_spectrum_power(power, biomass, per_log_size)
    # `power` on its own is still taken at face value, so plotCDF(power = 2)
    # keeps working. But a `power` that agrees with `biomass` only because the
    # spectrum is one in logarithmic size cannot be what the user meant.
    if (!is.null(biomass) && spectrum$per_log_size) {
        stop("You gave `power = ", power, "` together with `biomass = ",
             biomass, "`, which describes a density with respect to ",
             "logarithmic size. A cumulative distribution does not depend on ",
             "that choice, so use `biomass = ", biomass, "` on its own here. ",
             "It accumulates ", if (biomass) "biomass." else "numbers.")
    }
    spectrum
}

#' @rdname plotCDF
#' @usage NULL
#' @export
plotlyCDF <- function(object, species = NULL,
                      time_range, geometric_mean = FALSE,
                      wlim = c(NA, NA), llim = c(NA, NA),
                      ylim = c(NA, NA),
                      power = NULL, biomass = NULL, per_log_size = NULL,
                      total = FALSE, resource = FALSE,
                      background = TRUE,
                      highlight = NULL, normalise = TRUE,
                      log_x = TRUE, log_y = FALSE, log = NULL,
                      size_axis = c("w", "l"),
                      min_w = 0, max_w = Inf,
                      min_l = NULL, max_l = NULL, ...) {
    plotly_forward("plotCDF")
}

#' @rdname plotCDF2
#' @usage NULL
#' @export
plotlyCDF2 <- function(object1, object2, name1 = "First", name2 = "Second",
                       species = NULL,
                       wlim = c(NA, NA), llim = c(NA, NA),
                       ylim = c(NA, NA),
                       power = NULL, biomass = NULL, per_log_size = NULL,
                       total = FALSE, resource = FALSE,
                       background = TRUE,
                       highlight = NULL,
                       normalise = TRUE,
                       log_x = TRUE, log_y = FALSE,
                       log = NULL, size_axis = c("w", "l"),
                       min_w = 0, max_w = Inf,
                       min_l = NULL, max_l = NULL, ...) {
    plotly_forward("plotCDF2")
}
