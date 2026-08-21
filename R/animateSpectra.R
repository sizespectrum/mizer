#' Animate size-dependent quantities through time
#'
#' Creates an interactive plotly animation in which a play button steps through
#' time, drawing one line per species at each frame.
#'
#' The function dispatches on the class of `x`:
#'
#' * **`MizerSim`** — animates the community abundance spectra (number density
#'   or biomass density vs body size). Resource, background species, and a
#'   community total can be added via the `resource`, `background`, and `total`
#'   arguments. The `biomass` and `per_log_size` arguments choose the plotted
#'   quantity in the same way as in [plotSpectra()], and the `power` argument
#'   is available as the same alternative to them. Both axes are log10 by
#'   default and can each be switched to linear with `log_x = FALSE` or
#'   `log_y = FALSE`.
#'
#' * **`ArrayTimeBySpeciesBySize`** — animates any per-species, size-resolved
#'   quantity returned by a `MizerSim` accessor, such as [getFMort()],
#'   [getFeedingLevel()], or [getPredMort()]. Both axes are log10 by default
#'   and can each be switched to linear with `log_x = FALSE` or `log_y = FALSE`.
#'   Background species and a species total can be added via the `background`
#'   and `total` arguments.
#'
#' * **`ArrayTimeByResourceBySize`** — animates the size-resolved resource
#'   quantity returned by [NResource()] on a `MizerSim` object. There is only a
#'   single resource spectrum, so `species`, `total` and `background` do nothing
#'   and warn if set.
#'
#' Species linecolours and linetypes follow `params@linecolour` and
#' `params@linetype`.
#'
#' `animateSpectra()` is retained as a backward-compatible alias.
#'
#' @param x A `MizerSim`, `ArrayTimeBySpeciesBySize` or
#'   `ArrayTimeByResourceBySize` object.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted. Not used by the
#'   `ArrayTimeByResourceBySize` method, which warns if it is set.
#' @param tlim A numeric vector of length two providing lower and upper limits
#'   for the animated time window, e.g. `c(1997, 2007)`. Use `NA` to apply no
#'   limit at that end. Default is `c(NA, NA)`.
#' @param log_x If `TRUE` (default), use a log10 x-axis for body size.
#' @param log_y If `TRUE` (default), use a log10 y-axis.
#' @param log A character string specifying which axes to log-transform:
#'   `"x"`, `"y"`, `"xy"` or `""`. If supplied, this overrides `log_x`
#'   and `log_y`.
#' @param per_log_size Whether to animate a density per logarithmic size
#'   (`TRUE`) rather than per size (`FALSE`). Unlike `size_axis` this needs no
#'   weight-length relationship, so the `ArrayTimeByResourceBySize` method takes
#'   it too. The two kinds of `x` read the default `NULL` differently. For an
#'   array it means animating the density as it stands, and asking for it on an
#'   array that does not hold a density is an error. For a `MizerSim` it means
#'   `FALSE` unless `power` says otherwise, since there the argument chooses the
#'   plotted quantity together with `biomass` and `power`, in the same way as in
#'   [plotSpectra()].
#' @param size_axis Whether to plot size as weight (`"w"`, default) or length
#'   (`"l"`), using the allometric weight-length relationship of each species,
#'   or of the resource, see [resource_params()]. Number and biomass densities
#'   are transformed to match the chosen axis.
#' @param total A boolean value that determines whether the total over all
#'   selected species is plotted as an additional trace called `"Total"`.
#'   Default is `FALSE`. Not used by the `ArrayTimeByResourceBySize` method,
#'   which warns if it is set.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is `TRUE`. Not used by the `ArrayTimeByResourceBySize` method,
#'   which warns if it is set.
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the body-size (x) axis. Use `NA` to refer to the existing minimum or
#'   maximum.
#' @param llim A numeric vector of length two providing lower and upper limits
#'   for the length (x) axis when `size_axis = "l"`. Use `NA` to refer to the
#'   existing minimum or maximum.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the value (y) axis. Use `NA` to refer to the existing minimum or
#'   maximum. Limits are applied as Plotly axis ranges, so points outside the
#'   limits are clipped by the viewport rather than removed from the animation
#'   frames.
#' @param frame_duration Duration in milliseconds for which each saved frame is
#'   displayed. Default is 500.
#' @param transition_duration Duration in milliseconds of the interpolation
#'   between frames. Use `transition_duration = 0` to step directly from one
#'   saved frame to the next. Default is `frame_duration`.
#' @param easing The Plotly easing function to use when interpolating between
#'   frames. Default is `"linear"`. Available options are `"linear"`, `"quad"`,
#'   `"cubic"`, `"sin"`, `"exp"`, `"circle"`, `"elastic"`, `"back"`,
#'   `"bounce"`, and each of those with suffix `"-in"`, `"-out"`, or
#'   `"-in-out"` appended, for example `"cubic-in-out"`.
#' @param ... Further arguments used by only some of the methods:
#'
#'   **For `MizerSim` methods:**
#'   \describe{
#'     \item{`biomass`}{Whether to animate the biomass density (`TRUE`, the
#'       default) or the number density (`FALSE`).}
#'     \item{`power`}{The abundance is plotted as the number density times the
#'       weight raised to \code{power}. An alternative to `biomass` and
#'       `per_log_size`, with which it must agree if they are given as well.
#'       See [plotSpectra()] for details.}
#'     \item{`resource`}{A boolean value that determines whether resource is
#'       included. If `TRUE`, the resource spectrum is plotted as an additional
#'       trace called `"Resource"`. Default is `TRUE`.}
#'   }
#'
#' @return A plotly object with one animated line trace per plotted group. Use
#'   the play button or the slider to step through time.
#' @export
#' @family plotting functions
#' @examples
#' \donttest{
#' # Animate biomass density spectra, showing only sizes above 0.1 g
#' animate(NS_sim, power = 2, wlim = c(0.1, NA), tlim = c(1997, 2007))
#'
#' # Animate fishing mortality through time
#' animate(getFMort(NS_sim))
#'
#' # Animate feeding level for two species only
#' animate(getFeedingLevel(NS_sim), species = c("Cod", "Herring"))
#'
#' # Animate the resource spectrum
#' animate(NResource(NS_sim))
#' }
animate <- function(x, species = NULL, log_x = TRUE, log_y = TRUE,
                    log = NULL, wlim = c(NA, NA), llim = c(NA, NA),
                    ylim = c(NA, NA), tlim = c(NA, NA),
                    size_axis = c("w", "l"), per_log_size = NULL,
                    total = FALSE, background = TRUE, frame_duration = 500,
                    transition_duration = frame_duration,
                    easing = "linear", ...) {
    UseMethod("animate")
}

#' @rdname animate
#' @usage NULL
#' @export
animate.MizerSim <- function(x, species = NULL,
                              log_x = TRUE, log_y = TRUE,
                              log = NULL,
                              wlim = c(NA, NA), llim = c(NA, NA),
                              ylim = c(NA, NA),
                              tlim = c(NA, NA),
                              size_axis = c("w", "l"),
                              per_log_size = NULL,
                              total = FALSE,
                              background = TRUE,
                              frame_duration = 500,
                              transition_duration = frame_duration,
                              easing = "linear",
                              time_range = lifecycle::deprecated(),
                              power = NULL, biomass = NULL,
                              resource = TRUE, ...) {
    if (lifecycle::is_present(time_range)) {
        lifecycle::deprecate_warn("2.6.0", "animate(time_range)", "animate(tlim)")
        tlim <- c(min(time_range), max(time_range))
    }
    spectrum <- resolve_spectrum_power(power, biomass, per_log_size)
    power <- spectrum$power
    sim <- x
    size_axis <- plot_size_axis(size_axis)
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y
    assert_that(is.flag(total), is.flag(resource), is.flag(background),
                is.number(frame_duration), frame_duration >= 0,
                is.number(transition_duration), transition_duration >= 0,
                is.string(easing),
                length(wlim) == 2,
                length(llim) == 2,
                length(ylim) == 2)

    species <- valid_species_arg(sim, species)
    all_times <- as.numeric(dimnames(sim@n)$time)
    if (!is.na(tlim[1])) all_times <- all_times[all_times >= tlim[1]]
    if (!is.na(tlim[2])) all_times <- all_times[all_times <= tlim[2]]
    time_elements <- get_time_elements(sim, all_times)

    nf <- melt(sim@n[time_elements,
                     as.character(dimnames(sim@n)$sp) %in% species,
                               , drop = FALSE])
    names(nf)[names(nf) == "sp"] <- "Species"
    nf$legend_name <- as.character(nf$Species)

    # Add resource ----
    if (resource) {
        nf_pp <- melt(sim@n_pp[time_elements, , drop = FALSE])
        nf_pp$Species <- "Resource"
        nf_pp$legend_name <- "Resource"
        nf <- rbind(nf, nf_pp)
    }
    # Add background ----
    # Keep each background species as its own trace (avoids oscillation from
    # interleaved data points) but label them all as "Background" in the legend.
    if (background && any(sim@params@species_params$is_background)) {
        back_n <- sim@n[time_elements, sim@params@species_params$is_background, , drop = FALSE]
        nf_back <- melt(back_n)
        names(nf_back)[names(nf_back) == "sp"] <- "Species"
        nf_back$legend_name <- "Background"
        nf <- rbind(nf, nf_back)
    }
    # Add total ----
    # The contributors are assembled here but summed only after the size axis
    # has been converted, inside `animate_plotly()`. The total is the total of
    # everything the model holds: every species, whether or not it was selected
    # for display, and the resource, whether or not it is drawn — matching
    # `plotSpectra()`.
    total_dat <- NULL
    if (total) {
        total_dat <- melt(sim@n[time_elements, , , drop = FALSE])
        names(total_dat)[names(total_dat) == "sp"] <- "Species"
        nf_pp_total <- melt(sim@n_pp[time_elements, , drop = FALSE])
        nf_pp_total$Species <- "Resource"
        total_dat <- rbind(total_dat, nf_pp_total)
        total_dat$legend_name <- "Total"
    }

    y_label <- spectra_y_label(power, size_axis,
                               biomass = spectrum$biomass,
                               per_log_size = spectrum$per_log_size)
    # The animated spectrum is a density, so under second-order bin-averaging
    # we evaluate both the w^power weight and the plotted location at the
    # geometric bin centre w* = w sqrt(beta) (issue #383), matching
    # plotSpectra(). Default plots are unchanged.
    if (isTRUE(sim@params@second_order_w[["bin_average"]])) {
        beta <- sim@params@w_full[2] / sim@params@w_full[1]
        nf$w <- nf$w * sqrt(beta)
        if (!is.null(total_dat)) total_dat$w <- total_dat$w * sqrt(beta)
    }
    nf <- mutate(nf, value = value * w^power)
    if (!is.null(total_dat)) {
        total_dat <- mutate(total_dat, value = value * w^power)
    }

    animate_plotly(nf, sim@params, log_x, log_y, y_label, wlim, llim,
                   ylim,
                   size_axis = size_axis,
                   density_wrt = spectrum_density_wrt(spectrum$per_log_size),
                   total_dat = total_dat,
                   frame_duration = frame_duration,
                   transition_duration = transition_duration,
                   easing = easing)
}

# Build a plotly animation from a prepared long-format data frame.
# df must have columns: Species, legend_name, w, time, value.
# Traces are ordered by legend_name first (following params@linecolour), then
# by individual Species within each legend group — so background species always
# appear together and share a single legend entry.
animate_plotly <- function(df, params, log_x, log_y, y_label,
                           wlim = c(NA, NA), llim = c(NA, NA),
                           ylim = c(NA, NA),
                           size_axis = "w",
                           density_wrt = NA_character_,
                           per_log_size = NULL,
                           total_dat = NULL,
                           frame_duration = 500, transition_duration = 500,
                           easing = "linear") {
    size_axis <- plot_size_axis(size_axis)
    convert <- function(d) {
        convert_plot_density_axis(d, params, size_axis,
                                  density_wrt = density_wrt,
                                  per_log_size = per_log_size,
                                  value_col = "value")
    }
    df <- convert(df)
    x_var <- plot_size_x_var(size_axis)
    # The total is summed over the converted series, so that on a length axis
    # it is a sum at equal length rather than at equal weight.
    if (!is.null(total_dat)) {
        total_dat <- add_total_line(convert(total_dat), x_var, "value",
                                    by = "time")
        total_dat <- total_dat[total_dat$Species == "Total", ]
        total_dat$legend_name <- "Total"
        df <- rbind(df, total_dat[, names(df), drop = FALSE])
    }
    legend_name_order <- intersect(names(params@linecolour),
                                   unique(df$legend_name))
    sp_order <- unlist(lapply(legend_name_order, function(ln) {
        intersect(names(params@linecolour),
                  unique(df$Species[df$legend_name == ln]))
    }))
    p <- plotly::plot_ly()
    shown_legend_names <- character(0)
    for (sp in sp_order) {
        df_sp <- df[df$Species == sp, , drop = FALSE]
        ln <- unique(df_sp$legend_name)
        col <- params@linecolour[[ln]]
        showlegend <- !(ln %in% shown_legend_names)
        shown_legend_names <- c(shown_legend_names, ln)
        p <- plotly::add_lines(
            p,
            data = df_sp,
            x = stats::as.formula(paste0("~", x_var)), y = ~value,
            frame = ~time,
            name = ln,
            legendgroup = ln,
            line = list(color = col, simplify = FALSE),
            showlegend = showlegend
        )
    }
    p <- plotly::layout(p,
                        xaxis = plotly_axis(
                            df[[x_var]],
                            plot_size_xlim(wlim, size_axis, llim),
                            log_x, plot_size_xlab(size_axis)),
                        yaxis = plotly_axis(df$value, ylim, log_y, y_label),
                        legend = list(traceorder = "normal"))
    plotly::animation_opts(p, frame = frame_duration,
                           transition = transition_duration,
                           easing = easing)
}

plotly_axis <- function(values, limits, log_axis, title) {
    axis <- list(type = if (log_axis) "log" else "-",
                 exponentformat = "power",
                 title = title)
    range <- plotly_axis_range(values, limits, log_axis)
    if (!is.null(range)) axis$range <- range
    axis
}

plotly_axis_range <- function(values, limits, log_axis) {
    if (all(is.na(limits))) return(NULL)

    finite_values <- values[is.finite(values)]
    if (log_axis) finite_values <- finite_values[finite_values > 0]

    if (is.na(limits[1]) && length(finite_values) > 0) {
        limits[1] <- min(finite_values)
    }
    if (is.na(limits[2]) && length(finite_values) > 0) {
        limits[2] <- max(finite_values)
    }
    if (any(is.na(limits))) return(NULL)

    if (log_axis) {
        if (limits[1] <= 0 && length(finite_values) > 0) {
            limits[1] <- min(finite_values)
        }
        if (limits[2] <= 0 && length(finite_values) > 0) {
            limits[2] <- max(finite_values)
        }
        if (any(limits <= 0)) return(NULL)
        limits <- log10(limits)
    }

    if (!all(is.finite(limits))) return(NULL)
    if (limits[1] == limits[2]) {
        delta <- if (log_axis) 0.5 else max(abs(limits[1]) * 0.05, 1)
        limits <- limits + c(-1, 1) * delta
    }
    limits
}

#' @rdname animate
#' @param sim A `MizerSim` object.
#' @export
animateSpectra <- function(sim, ...) animate(sim, ...)
