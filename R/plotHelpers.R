# Shared size-axis, density-conversion, proportion and total helpers.

#' Validate the size-axis argument
#'
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @return The matched size-axis string, either `"w"` or `"l"`.
#' @keywords internal
plot_size_axis <- function(size_axis = "w") {
    match.arg(size_axis, c("w", "l"))
}

#' Name of the x-variable for a given size axis
#'
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @return `"l"` for a length axis, otherwise `"w"`.
#' @keywords internal
plot_size_x_var <- function(size_axis) {
    if (identical(plot_size_axis(size_axis), "l")) "l" else "w"
}

#' Axis label for a given size axis
#'
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @return The axis label: `"Length [cm]"` for a length axis, otherwise
#'   `"Size [g]"`.
#' @keywords internal
plot_size_xlab <- function(size_axis) {
    if (identical(plot_size_axis(size_axis), "l")) "Length [cm]" else "Size [g]"
}

#' Choose the x-axis limits for a given size axis
#'
#' @param wlim Numeric vector of length two giving the weight limits.
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @param llim Numeric vector of length two giving the length limits.
#' @return `llim` for a length axis, otherwise `wlim`.
#' @keywords internal
plot_size_xlim <- function(wlim, size_axis, llim = c(NA, NA)) {
    if (identical(plot_size_axis(size_axis), "l")) llim else wlim
}

#' Assemble the tooltip variables for a size-axis plot
#'
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @param before,after Optional character vectors of variable names to place
#'   before and after the size variable in the tooltip.
#' @return A character vector of tooltip variable names.
#' @keywords internal
plot_size_tooltip <- function(size_axis, before = NULL, after = NULL) {
    c(before, plot_size_x_var(size_axis), after)
}

#' Y-axis limits for a plot of a proportion
#'
#' A proportion is easiest to read against the whole of the interval from 0 to
#' 1, so that is the range a plot of one shows by default. The range is only
#' ever *widened* to include the data, never narrowed to the interval: a
#' critical feeding level or a resource level above 1 is a real feature of the
#' model and must stay visible.
#'
#' Only the ends of `ylim` that the caller left as `NA` are filled in, so an
#' explicit limit always wins. A logarithmic axis is left alone, having no place
#' for the 0.
#'
#' @param ylim Numeric vector of length two, the limits the caller asked for.
#' @param log_y Whether the y axis is logarithmic.
#' @param values The values being plotted.
#' @return A numeric vector of length two.
#' @keywords internal
proportion_ylim <- function(ylim, log_y, values) {
    if (isTRUE(log_y)) return(ylim)
    values <- values[is.finite(values)]
    full <- range(c(0, 1, values))
    if (is.na(ylim[1])) ylim[1] <- full[1]
    if (is.na(ylim[2])) ylim[2] <- full[2]
    ylim
}

#' Whether a plot's y axis is logarithmic
#'
#' Read from the plot's own y scale, so that data being *added* to an existing
#' plot can be filtered the way that plot was. A plot with no explicit y scale
#' has ggplot2's default, which is linear.
#'
#' @param plot A ggplot2 object.
#' @return `TRUE` if the plot's y axis uses a log10 transformation.
#' @keywords internal
plot_y_is_log <- function(plot) {
    for (scale in plot$scales$scales) {
        if ("y" %in% scale$aesthetics) {
            return(identical(scale$trans$name, "log-10"))
        }
    }
    FALSE
}

#' Y-axis limits an array's type calls for
#'
#' Only a `"proportion"` has an opinion. A `"density"` is handled where the size
#' axis is converted, and a `"value"` needs nothing.
#'
#' @param x A mizer array object.
#' @inheritParams proportion_ylim
#' @return A numeric vector of length two.
#' @keywords internal
array_ylim <- function(x, ylim, log_y, values) {
    if (identical(array_type(x), "proportion")) {
        return(proportion_ylim(ylim, log_y, values))
    }
    ylim
}

#' The logarithmic y axis an array's type calls for
#'
#' A proportion belongs on a linear axis, so a plot of one turns `log_y` off
#' unless the caller asked for a particular axis. Called before
#' [parsePlotLog()], which is why it also has to check `log`.
#'
#' @param x A mizer array object.
#' @param log_y The `log_y` argument of the plot method.
#' @param log The `log` argument of the plot method.
#' @param given Whether the caller supplied `log_y` (i.e. `!missing(log_y)`).
#' @return The `log_y` to use.
#' @keywords internal
array_log_y <- function(x, log_y, log, given) {
    if (!given && is.null(log) && identical(array_type(x), "proportion")) {
        return(FALSE)
    }
    log_y
}

#' Add a total line to plotting data by summing over its series
#'
#' The total has to be formed *after* the size coordinate has been converted,
#' not before. On a weight axis every series shares the model's weight grid, so
#' summing at equal weight and summing at equal position are the same thing. On
#' a length axis they are not: each species, and the resource, converts weight
#' to length with its own allometric relationship, so at a given length the
#' series sit at different weights and their grids no longer coincide. The sum
#' that means something there is the sum at equal *length* — the number of
#' organisms per unit length, whatever they are — which is what this computes.
#'
#' Each series is interpolated onto the sorted union of all the size
#' coordinates, linearly in the logarithm of size, since the grid is
#' logarithmic. A density contributes nothing outside its own range; a
#' cumulative series retains its final value beyond its range. When the
#' series already share a grid — always on a weight axis, and on a length axis
#' whenever the weight-length parameters agree — the union is that grid and the
#' interpolation reproduces the values exactly, so nothing is approximated in
#' the cases where nothing needs to be.
#'
#' @param plot_dat A data frame of plotting data with a size column, a value
#'   column and a `Species` column.
#' @param x_var Name of the size column. Defaults to the first column.
#' @param value_col Name or index of the value column. Defaults to the second.
#' @param by Names of further columns identifying separate plots, such as the
#'   time of an animation frame or the model of a comparison. A total is formed
#'   within each of their combinations.
#' @param cumulative Whether completed cumulative curves retain their final
#'   value beyond their last size, rather than contributing zero there.
#' @return `plot_dat` with the total appended as a series named `"Total"`.
#' @keywords internal
add_total_line <- function(plot_dat, x_var = names(plot_dat)[[1]],
                           value_col = 2, by = NULL, cumulative = FALSE) {
    if (nrow(plot_dat) == 0) return(plot_dat)
    if (is.numeric(value_col)) value_col <- names(plot_dat)[[value_col]]
    if (length(by) > 0) {
        groups <- split(plot_dat, plot_dat[, by, drop = FALSE], drop = TRUE)
        out <- lapply(groups, add_total_line, x_var = x_var,
                      value_col = value_col, cumulative = cumulative)
        out <- do.call(rbind, unname(out))
        rownames(out) <- NULL
        return(out)
    }
    x_all <- sort(unique(plot_dat[[x_var]]))
    total <- numeric(length(x_all))
    for (series in split(plot_dat, plot_dat$Species)) {
        series <- series[order(series[[x_var]]), , drop = FALSE]
        x <- series[[x_var]]
        y <- series[[value_col]]
        if (length(x) < 2) {
            # `approx()` needs two points. A single density value contributes
            # only at that coordinate; a cumulative value persists beyond it.
            if (cumulative) {
                total[x_all >= x] <- total[x_all >= x] + y
            } else {
                hit <- match(x, x_all)
                total[hit] <- total[hit] + y
            }
            next
        }
        contribution <- stats::approx(log(x), y, xout = log(x_all),
                                      yleft = if (cumulative) 0 else NA_real_,
                                      yright = if (cumulative) tail(y, 1) else
                                          NA_real_)$y
        contribution[is.na(contribution)] <- 0
        total <- total + contribution
    }
    total_dat <- plot_dat[rep(1, length(x_all)), , drop = FALSE]
    total_dat[[x_var]] <- x_all
    total_dat[[value_col]] <- total
    total_dat$Species <- "Total"
    if ("Legend" %in% names(total_dat)) total_dat$Legend <- "Total"
    rownames(total_dat) <- NULL
    rbind(plot_dat, total_dat)
}

#' Convert the contributors to a total and append the total line
#'
#' Wraps the two steps that every plot of an array with `total = TRUE` needs:
#' the contributors are put on the same size axis as the data they will join,
#' and then summed there, see [add_total_line()].
#'
#' @param plot_dat The converted plotting data to append the total to.
#' @param total_dat The unconverted contributors, see [total_contributors()].
#' @param params A MizerParams object.
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @param x The mizer array the data came from, which says whether the values
#'   are a density.
#' @param per_log_size Whether to express a density per logarithmic size.
#' @return `plot_dat` with the total appended as a series named `"Total"`.
#' @keywords internal
append_total_line <- function(plot_dat, total_dat, params, size_axis, x,
                              per_log_size = NULL) {
    total_dat <- convert_plot_density_axis(total_dat, params, size_axis,
                                           density_wrt = array_density_wrt(x),
                                           per_log_size = per_log_size)
    total_dat <- add_total_line(total_dat)
    total_dat <- total_dat[total_dat$Species == "Total", , drop = FALSE]
    rbind(plot_dat, total_dat[, names(plot_dat), drop = FALSE])
}

#' The weight-length parameters to plot each row of plotting data with
#'
#' A length axis needs an allometric weight-length relationship for every line
#' on the plot. The species take theirs from their species parameters and the
#' resource takes its from [resource_params()], where it defaults to the
#' equivalent spherical diameter (see [resource_length_params()]). Anything else
#' — the "Total" row, for instance — has none, and is reported as `NA` so that
#' the caller can leave it out.
#'
#' @param species A vector of the species names in the plotting data.
#' @param params A MizerParams object providing the weight-length parameters.
#' @return A data frame with columns `a` and `b`, one row for each element of
#'   `species`, holding `NA` where no relationship is known.
#' @keywords internal
plot_length_params <- function(species, params) {
    sp <- params@species_params
    if (!all(c("a", "b") %in% names(sp))) {
        stop("The species parameter data frame must contain columns ",
             "'a' and 'b'.")
    }
    species <- as.character(species)
    idx <- match(species, as.character(sp$species))
    out <- data.frame(a = sp[["a"]][idx], b = sp[["b"]][idx])
    is_resource <- species == "Resource"
    if (any(is_resource)) {
        rp <- resource_length_params(params)
        out$a[is_resource] <- rp$a
        out$b[is_resource] <- rp$b
    }
    out
}

#' Convert plotting data from weight to length
#'
#' When `size_axis = "l"`, adds a length column `l` computed from the weight
#' column `w` using the weight-length relationship of each line, see
#' [plot_length_params()]. Rows with no such relationship are dropped. For
#' `size_axis = "w"` the data is returned unchanged.
#'
#' @param plot_dat A data frame of plotting data with a `w` column and a species
#'   column.
#' @param params A MizerParams object providing the weight-length parameters.
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @param species_col Name of the column identifying the species. Default is
#'   `"Species"`.
#' @param drop_w If `TRUE` (default), the `w` column is dropped once `l` has been
#'   computed.
#' @return The plotting data, with a length column `l` added (and `w` optionally
#'   dropped) when `size_axis = "l"`.
#' @keywords internal
convert_plot_size_axis <- function(plot_dat, params, size_axis,
                                   species_col = "Species",
                                   drop_w = TRUE) {
    size_axis <- plot_size_axis(size_axis)
    if (identical(size_axis, "w")) {
        return(plot_dat)
    }
    if (!is(params, "MizerParams")) {
        stop("Length-axis plots need a MizerParams object.")
    }
    if (!"w" %in% names(plot_dat)) {
        stop("Length-axis plots need plotting data with a `w` column.")
    }
    if (!species_col %in% names(plot_dat)) {
        stop("Length-axis plots need plotting data with a `", species_col,
             "` column.")
    }

    ab <- plot_length_params(plot_dat[[species_col]], params)
    known <- !is.na(ab$a) & !is.na(ab$b)
    plot_dat <- plot_dat[known, , drop = FALSE]
    ab <- ab[known, , drop = FALSE]
    if (nrow(plot_dat) == 0) {
        if (!drop_w) {
            plot_dat$l <- numeric(0)
            return(plot_dat[, c("l", "w",
                                setdiff(names(plot_dat), c("l", "w"))),
                            drop = FALSE])
        }
        return(plot_dat[, setdiff(names(plot_dat), "w"), drop = FALSE])
    }
    plot_dat$l <- (plot_dat$w / ab$a)^(1 / ab$b)
    if (drop_w) {
        return(plot_dat[, c("l", setdiff(names(plot_dat), c("l", "w"))),
                        drop = FALSE])
    }
    plot_dat[, c("l", "w", setdiff(names(plot_dat), c("l", "w"))),
             drop = FALSE]
}

#' Density measures a spectrum can be expressed in
#'
#' A size spectrum is a density, and a density only has a meaning together with
#' the variable it is a density with respect to. That variable is one of
#' \describe{
#'   \item{`"w"`}{a density with respect to weight, e.g. numbers per gram.}
#'   \item{`"log_w"`}{a density with respect to logarithmic weight, e.g.
#'     numbers per log weight interval.}
#'   \item{`"l"`}{a density with respect to length, e.g. numbers per cm.}
#'   \item{`"log_l"`}{a density with respect to logarithmic length.}
#'   \item{`NA`}{not a density, e.g. a rate or a dimensionless quantity. Such
#'     values are left alone when the size axis changes.}
#' }
#'
#' Mizer arrays are indexed by the model's weight grid, so an array that holds
#' a density (`type = "density"`, see [array_types]) always holds one with
#' respect to weight. The other measures arise for quantities the spectrum plots
#' compute on the fly: `plotSpectra(per_log_size = TRUE)` shows a density with
#' respect to logarithmic weight, and either can be restated per unit length by
#' `size_axis = "l"`.
#'
#' @format A character vector of the four density measures.
#' @keywords internal
density_measures <- c("w", "log_w", "l", "log_l")

#' Validate a density measure
#'
#' @param density_wrt A density measure, see [density_measures]. `NULL` and `NA`
#'   both stand for "not a density".
#' @return The validated measure, or `NA_character_` when the values are not a
#'   density.
#' @keywords internal
validate_density_wrt <- function(density_wrt) {
    if (is.null(density_wrt)) return(NA_character_)
    if (length(density_wrt) == 1 && all(is.na(density_wrt))) {
        return(NA_character_)
    }
    if (!is.character(density_wrt) || length(density_wrt) != 1 ||
            !density_wrt %in% density_measures) {
        stop("`density_wrt` must be NA or one of ",
             paste0("\"", density_measures, "\"", collapse = ", "), ".")
    }
    density_wrt
}

#' The density measure a plot calls for
#'
#' A density is expressed with respect to two independent choices: the size
#' variable, which follows `size_axis`, and whether it is per size or per
#' logarithmic size, which follows `per_log_size`. Plotting against a length
#' axis therefore turns a density with respect to weight into one with respect
#' to length, and a density with respect to logarithmic weight into one with
#' respect to logarithmic length.
#'
#' @param density_wrt The measure the values are a density with respect to, see
#'   [density_measures].
#' @param size_axis Either `"w"` (weight) or `"l"` (length).
#' @param per_log_size Whether to express the values per logarithmic size.
#'   `NULL` (the default) keeps whichever the values already are.
#' @return The density measure to express the values in, or `NA_character_` if
#'   the values are not a density.
#' @keywords internal
density_target_measure <- function(density_wrt, size_axis,
                                   per_log_size = NULL) {
    density_wrt <- validate_density_wrt(density_wrt)
    if (is.na(density_wrt)) return(NA_character_)
    per_log_size <- per_log_size %||% startsWith(density_wrt, "log_")
    paste0(if (isTRUE(per_log_size)) "log_" else "", plot_size_axis(size_axis))
}

#' Factor relating a density measure to a density with respect to weight
#'
#' Writing \eqn{N_w} for the density with respect to weight, the density with
#' respect to measure \eqn{m} is \eqn{N_w} times the factor returned here. With
#' the allometric weight-length relationship \eqn{w = a l^b} these factors are
#' \eqn{1} for `"w"`, \eqn{w} for `"log_w"`, \eqn{dw/dl = b w / l} for `"l"`
#' and \eqn{l\,dw/dl = b w} for `"log_l"`.
#'
#' @param measure One of [density_measures].
#' @param w,l,b Numeric vectors of the same length giving the weight, the
#'   corresponding length, and the exponent of the weight-length relationship.
#' @return A numeric vector of factors.
#' @keywords internal
density_measure_weight <- function(measure, w, l, b) {
    switch(measure,
           w = rep(1, length(w)),
           log_w = w,
           l = b * w / l,
           log_l = b * w,
           stop("Unknown density measure `", measure, "`."))
}

#' Jacobian converting between two density measures
#'
#' @param from,to Density measures, see [density_measures].
#' @inheritParams density_measure_weight
#' @return A numeric vector by which to multiply a density with respect to
#'   `from` to obtain the density with respect to `to`.
#' @keywords internal
density_measure_jacobian <- function(from, to, w, l, b) {
    if (identical(from, to)) return(rep(1, length(w)))
    density_measure_weight(to, w, l, b) / density_measure_weight(from, w, l, b)
}

#' The size unit appearing in the units of a density
#'
#' @param measure One of [density_measures].
#' @return `"g"` or `"cm"`, or `NA_character_` for a density with respect to a
#'   logarithmic size, whose units carry no size unit.
#' @keywords internal
density_size_unit <- function(measure) {
    switch(measure, w = "g", l = "cm", NA_character_)
}

#' Restate the units of a density in a different density measure
#'
#' The per-size factor is found in either of the two spellings mizer uses for
#' it, `1/g` and `g^-1`, and is then swapped for the size unit of the target
#' measure. A density per *logarithmic* size carries no size unit at all — a
#' number per log weight interval is just a number — so converting to one
#' removes the factor instead: `1/g` becomes dimensionless and `g^-1/year`
#' becomes `1/year`.
#'
#' Units that state no per-size factor are returned unchanged, since there is
#' then nothing to identify as the size unit.
#'
#' @param units The units of the values, possibly `NULL`.
#' @param from,to Density measures, see [density_measures].
#' @return The units expressed in the `to` measure.
#' @keywords internal
convert_density_units <- function(units, from, to) {
    if (is.na(from) || is.na(to) || identical(from, to)) return(units)
    to_sym <- density_size_unit(to)
    if (is.null(units) || !nzchar(units)) {
        # Nothing was declared, so fall back to the canonical density units
        return(if (is.na(to_sym)) units else paste0("1/", to_sym))
    }
    from_sym <- density_size_unit(from)
    if (is.na(from_sym)) return(units)
    # Replacing the per-size factor by "1" removes it; a units string that is
    # then nothing but "1" is better shown as no units at all.
    replacement <- if (is.na(to_sym)) "1" else paste0("1/", to_sym)
    units <- sub(paste0("\\b1/", from_sym, "\\b"), replacement, units)
    replacement <- if (is.na(to_sym)) "1" else paste0(to_sym, "^-1")
    units <- sub(paste0("\\b", from_sym, "\\^-1"), replacement, units)
    if (identical(units, "1")) "" else units
}

#' Express plotting data on the requested size axis
#'
#' Converts the size coordinate of the plotting data to the requested axis and,
#' when the values are a density, multiplies them by the Jacobian that restates
#' them in the density measure that axis calls for (see
#' [density_target_measure()]). Values that are not a density are left alone.
#'
#' Anything involving a length is a per-species quantity, because the
#' weight-length relationship is, so rows whose species is not one of the
#' model's species — the "Total" row, for instance — are dropped when a length
#' is needed. Going from a density per size to one per logarithmic size needs no
#' length, and keeps those rows.
#'
#' @inheritParams convert_plot_size_axis
#' @param density_wrt The measure the values are a density with respect to, see
#'   [density_measures]. `NA` (the default) means the values are not a density.
#' @param per_log_size Whether to express the values per logarithmic size.
#'   `NULL` (the default) keeps whichever the values already are.
#' @param value_col Name or index of the value column. Defaults to the second
#'   column.
#' @return The plotting data with its size coordinate, and where called for its
#'   values, expressed for the requested axis. The size coordinate is the first
#'   column.
#' @keywords internal
convert_plot_density_axis <- function(plot_dat, params, size_axis,
                                      density_wrt = NA_character_,
                                      per_log_size = NULL,
                                      species_col = "Species",
                                      value_col = 2) {
    size_axis <- plot_size_axis(size_axis)
    density_wrt <- validate_density_wrt(density_wrt)
    target <- density_target_measure(density_wrt, size_axis, per_log_size)
    needs_jacobian <- !is.na(target) && !identical(target, density_wrt)
    if (identical(size_axis, "w") && !needs_jacobian) {
        return(plot_dat)
    }
    if (is.numeric(value_col)) {
        value_col <- names(plot_dat)[[value_col]]
    }
    # A length axis needs the lengths, and so does any Jacobian that involves
    # one. The weights are kept because the Jacobians are written in terms of
    # both.
    needs_length <- identical(size_axis, "l") ||
        isTRUE(density_wrt %in% c("l", "log_l"))
    if (needs_length) {
        plot_dat <- convert_plot_size_axis(plot_dat, params, "l",
                                           species_col = species_col,
                                           drop_w = FALSE)
    }
    if (needs_jacobian && nrow(plot_dat) > 0) {
        if (needs_length) {
            l <- plot_dat$l
            b <- plot_length_params(plot_dat[[species_col]], params)$b
        } else {
            # Unused by a Jacobian between two weight-based measures. Passing
            # NA rather than a placeholder keeps a mistake here loud.
            l <- b <- rep(NA_real_, nrow(plot_dat))
        }
        plot_dat[[value_col]] <- plot_dat[[value_col]] *
            density_measure_jacobian(density_wrt, target, plot_dat$w, l, b)
    }
    x_var <- plot_size_x_var(size_axis)
    plot_dat[, c(x_var, setdiff(names(plot_dat), c("l", "w"))), drop = FALSE]
}
