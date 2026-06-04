# ArrayTimeBySpecies S3 class for time x species arrays
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' S3 class for time x species arrays
#'
#' Some functions in mizer return two-dimensional arrays (time x species)
#' holding quantities like biomass, abundance, or yield rate through time.
#' The `ArrayTimeBySpecies` class wraps these arrays to provide convenient
#' `print()`, `summary()`, `plot()`, and `as.data.frame()` methods.
#'
#' An `ArrayTimeBySpecies` object behaves just like a regular matrix for
#' arithmetic operations and subsetting. It carries these lightweight attributes:
#' \itemize{
#'   \item `value_name` – a human-readable name for the value
#'       (e.g. "Biomass").
#'   \item `units` – the units of the value (e.g. "g", "g/year").
#'   \item `params` – the `MizerParams` object that created the values.
#' }
#'
#' @param x A matrix (time x species).
#' @param value_name A string giving the human-readable name for the value.
#' @param units A string giving the units (e.g. "g", "g/year").
#' @param params A `MizerParams` object holding the model that created the
#'   values.
#'
#' @return An `ArrayTimeBySpecies` object (inherits from `matrix` and `array`).
#'
#' @seealso [print()], [summary()], [as.data.frame()], [plot()]
#' @export
#' @examples
#' \donttest{
#' bio <- getBiomass(NS_sim)
#' is.ArrayTimeBySpecies(bio)
#' summary(bio)
#' }
ArrayTimeBySpecies <- function(x, value_name = NULL, units = NULL,
                               params = NULL) {
    if (!is.matrix(x)) {
        stop("`x` must be a matrix.")
    }
    structure(x,
        class = c("ArrayTimeBySpecies", "matrix", "array"),
        value_name = value_name,
        units = units,
        params = params
    )
}

#' Test if an object is a ArrayTimeBySpecies
#'
#' @param x An object to test.
#' @return `TRUE` if `x` is an `ArrayTimeBySpecies` object, `FALSE` otherwise.
#' @export
#' @examples
#' is.ArrayTimeBySpecies(getBiomass(NS_sim))
#' is.ArrayTimeBySpecies(matrix(1:4, nrow = 2))
is.ArrayTimeBySpecies <- function(x) {
    inherits(x, "ArrayTimeBySpecies")
}

#' @export
print.ArrayTimeBySpecies <- function(x, ...) {
    value_name <- attr(x, "value_name") %||% "ArrayTimeBySpecies"
    units_str <- attr(x, "units")
    dims <- dim(x)
    header <- paste0(value_name, " (", dims[1], " times x ", dims[2],
                     " species)")
    if (!is.null(units_str) && nzchar(units_str)) {
        header <- paste0(header, " [", units_str, "]")
    }
    cat(header, "\n")
    sp_names <- colnames(x)
    if (!is.null(sp_names)) {
        vals <- apply(unclass(x), 2, function(col) {
            col <- col[is.finite(col)]
            if (length(col) == 0) return("all NA/Inf")
            paste0("min=", signif(min(col), 3),
                   " mean=", signif(mean(col), 3),
                   " max=", signif(max(col), 3))
        })
        for (i in seq_along(sp_names)) {
            cat("  ", sp_names[i], ": ", vals[i], "\n", sep = "")
        }
    }
    invisible(x)
}

#' @export
summary.ArrayTimeBySpecies <- function(object, ...) {
    value_name <- attr(object, "value_name") %||% "ArrayTimeBySpecies"
    units_str <- attr(object, "units")
    sp_names <- colnames(object)
    mat <- unclass(object)

    df <- data.frame(
        Species = sp_names,
        Min = apply(mat, 2, min, na.rm = TRUE),
        Mean = apply(mat, 2, mean, na.rm = TRUE),
        Max = apply(mat, 2, max, na.rm = TRUE),
        row.names = NULL,
        stringsAsFactors = FALSE
    )

    result <- list(
        value_name = value_name,
        units = units_str,
        dims = dim(object),
        per_species = df
    )
    class(result) <- "summary.ArrayTimeBySpecies"
    result
}

#' @export
print.summary.ArrayTimeBySpecies <- function(x, ...) {
    header <- x$value_name
    if (!is.null(x$units) && nzchar(x$units)) {
        header <- paste0(header, " [", x$units, "]")
    }
    cat(header, "\n")
    cat(x$dims[1], "times x", x$dims[2], "species\n\n")
    print(x$per_species, row.names = FALSE)
    invisible(x)
}

#' @rdname plot
#' @usage NULL
#'
#' @param start_time The first time to be plotted. Default (`NULL`) is the
#'   beginning of the time series. Only applies to `ArrayTimeBySpecies`.
#' @param end_time The last time to be plotted. Default (`NULL`) is the end of
#'   the time series. Only applies to `ArrayTimeBySpecies`.
#' @export
#' @examples
#' \donttest{
#' plot(getBiomass(NS_sim))
#' plot(getBiomass(NS_sim), species = c("Cod", "Herring"), total = TRUE)
#' plot(getYield(NS_sim), species = c("Cod", "Herring"))
#' }
plot.ArrayTimeBySpecies <- function(x, species = NULL,
                                    start_time = NULL, end_time = NULL,
                                    y_ticks = 6, ylim = c(NA, NA),
                                    total = FALSE, background = TRUE,
                                    highlight = NULL, log_x = FALSE,
                                    log_y = TRUE, log = NULL,
                                    return_data = FALSE,
                                    ...) {
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    value_name <- attr(x, "value_name") %||% "Value"
    units_str <- attr(x, "units")
    params <- attr(x, "params")

    # Construct y_label
    y_label <- value_name
    if (!is.null(units_str) && nzchar(units_str)) {
        y_label <- paste0(value_name, " [", units_str, "]")
    }
    plot_dat <- prepare_ArrayTimeBySpecies_plot_data(
        x, species = species, start_time = start_time, end_time = end_time,
        ylim = ylim, total = total, background = background)

    if (return_data) return(plot_dat)

    plotDataFrame(plot_dat, params, xlab = "Year", ylab = y_label,
                  xtrans = if (log_x) "log10" else "identity",
                  ytrans = if (log_y) "log10" else "identity",
                  ylim = ylim, y_ticks = y_ticks, highlight = highlight,
                  legend_var = "Legend")
}

#' @rdname addPlot
#' @usage NULL
#' @inheritParams plot
#' @param colour Optional fixed colour for the added lines. If `NULL`, the
#'   species colours from the existing plot are used.
#' @param linetype Optional fixed line type for the added lines. If `NULL`, the
#'   species line types from the existing plot are used.
#' @param linewidth Width of the added lines.
#' @param alpha Transparency of the added lines.
#' @export
addPlot.ArrayTimeBySpecies <- function(plot, x, species = NULL,
                                       start_time = NULL,
                                       end_time = NULL,
                                       ylim = c(NA, NA),
                                       total = FALSE,
                                       background = TRUE,
                                       colour = NULL,
                                       linetype = "dashed",
                                       linewidth = 0.8,
                                       alpha = 1,
                                       ...) {
    if (!inherits(plot, "ggplot")) {
        stop("The `plot` argument must be a ggplot object.")
    }
    assert_that(is.number(linewidth),
                is.number(alpha),
                alpha >= 0,
                alpha <= 1)

    plot <- deep_copy(plot)
    plot_dat <- prepare_ArrayTimeBySpecies_plot_data(
        x, species = species, start_time = start_time, end_time = end_time,
        ylim = ylim, total = total, background = background)
    y_var <- names(plot_dat)[[2]]
    check_addPlot_compatible(plot, x_var = "Year", y_var = y_var,
                             units = attr(x, "units"))

    mapping <- aes(x = .data[["Year"]], y = .data[[y_var]],
                   group = .data[["Species"]])
    if (is.null(colour)) {
        mapping$colour <- rlang::quo(.data[["Legend"]])
    }
    if (is.null(linetype)) {
        mapping$linetype <- rlang::quo(.data[["Legend"]])
    }

    layer_args <- list(
        data = plot_dat,
        mapping = mapping,
        linewidth = linewidth,
        alpha = alpha,
        inherit.aes = FALSE
    )
    if (!is.null(colour)) {
        layer_args$colour <- colour
    }
    if (!is.null(linetype)) {
        layer_args$linetype <- linetype
    }

    plot + do.call(geom_line, layer_args)
}

#' @rdname plot
#' @usage NULL
#' @export
plot2.ArrayTimeBySpecies <- function(x, y, name1 = "First", name2 = "Second",
                                     species = NULL,
                                     start_time = NULL, end_time = NULL,
                                     y_ticks = 6, ylim = c(NA, NA),
                                     total = FALSE, background = TRUE,
                                     log_x = FALSE, log_y = TRUE,
                                     log = NULL, ...) {
    check_plot2_compatible(x, y, "ArrayTimeBySpecies")
    compare_array_metadata(x, y)
    log_axes <- parsePlotLog(log, log_x = log_x, log_y = log_y)
    log_x <- log_axes$log_x
    log_y <- log_axes$log_y

    params <- attr(x, "params")
    y_label <- array_y_label(x, default = "Value")
    plot_dat1 <- prepare_ArrayTimeBySpecies_plot_data(
        x, species = species, start_time = start_time, end_time = end_time,
        ylim = ylim, total = total, background = background)
    plot_dat2 <- prepare_ArrayTimeBySpecies_plot_data(
        y, species = species, start_time = start_time, end_time = end_time,
        ylim = ylim, total = total, background = background)

    plotComparisonDataFrame(plot_dat1, plot_dat2, params,
                            name1 = name1, name2 = name2,
                            xlab = "Year", ylab = y_label,
                            xtrans = if (log_x) "log10" else "identity",
                            ytrans = if (log_y) "log10" else "identity",
                            ylim = ylim, y_ticks = y_ticks,
                            legend_var = "Legend")
}

#' @rdname plot
#' @usage NULL
#' @export
plotRelative.ArrayTimeBySpecies <- function(x, y, species = NULL,
                                            start_time = NULL,
                                            end_time = NULL,
                                            ylim = c(NA, NA),
                                            total = FALSE,
                                            background = TRUE,
                                            log_x = FALSE, ...) {
    check_plot2_compatible(x, y, "ArrayTimeBySpecies")
    compare_array_metadata(x, y)
    params <- attr(x, "params")
    plot_dat1 <- prepare_ArrayTimeBySpecies_plot_data(
        x, species = species, start_time = start_time, end_time = end_time,
        total = total, background = background)
    plot_dat2 <- prepare_ArrayTimeBySpecies_plot_data(
        y, species = species, start_time = start_time, end_time = end_time,
        total = total, background = background)

    plotRelativeDataFrame(plot_dat1, plot_dat2, params,
                          xlab = "Year",
                          xtrans = if (log_x) "log10" else "identity",
                          ylim = ylim, legend_var = "Legend")
}

prepare_ArrayTimeBySpecies_plot_data <- function(x, species = NULL,
                                                 start_time = NULL,
                                                 end_time = NULL,
                                                 ylim = c(NA, NA),
                                                 total = FALSE,
                                                 background = TRUE) {
    value_name <- attr(x, "value_name") %||% "Value"
    params <- attr(x, "params")

    # Filter to time range
    t <- as.numeric(rownames(x))
    if (any(is.na(t))) t <- seq_len(nrow(x))
    if (!is.null(start_time) && !is.null(end_time) && start_time >= end_time) {
        stop("start_time must be less than end_time")
    }
    if (!is.null(start_time)) {
        x <- x[t >= start_time, , drop = FALSE]
        t <- t[t >= start_time]
    }
    if (!is.null(end_time)) {
        x <- x[t <= end_time, , drop = FALSE]
    }

    all_species <- colnames(x)
    if (!is.null(species)) {
        species <- intersect(species, all_species)
        if (length(species) == 0) {
            stop("None of the selected species are in the array.")
        }
    } else {
        species <- all_species
    }

    bm <- unclass(x)

    # Add total if requested (before species selection so total includes all)
    if (total) {
        bm <- cbind(bm, Total = rowSums(bm))
    }

    bm <- reshape2::melt(bm)

    # Implement ylim and a minimal cutoff; bring columns in desired order
    min_value <- 1e-20
    bm <- bm[bm$value >= min_value &
                 (is.na(ylim[1]) | bm$value >= ylim[1]) &
                 (is.na(ylim[2]) | bm$value <= ylim[2]), c(1, 3, 2)]
    names(bm) <- c("Year", value_name, "Species")

    # Select species
    plot_dat <- bm[bm$Species %in% c("Total", species), ]
    plot_dat$Legend <- plot_dat$Species

    if (background && !is.null(params) && any(params@species_params$is_background)) {
        # Add background species in light grey
        bkgrd_sp <- params@species_params$species[params@species_params$is_background]
        if (length(bkgrd_sp) > 0) {
            bm_bkgrd <- bm[bm$Species %in% bkgrd_sp, ]
            bm_bkgrd$Legend <- "Background"
            plot_dat <- rbind(plot_dat, bm_bkgrd)
        }
    }

    plot_dat
}

#' @rdname plot
#' @usage NULL
#' @exportS3Method plotly::ggplotly
#' @examples
#' \donttest{
#' ggplotly(getBiomass(NS_sim))
#' }
ggplotly.ArrayTimeBySpecies <- function(p = ggplot2::last_plot(),
                                        width = NULL, height = NULL,
                                        tooltip = "all",
                                        dynamicTicks = FALSE,
                                        layerData = 1,
                                        originalData = TRUE,
                                        source = "A", ...) {
    plotly::ggplotly(plot(p, ...), width = width, height = height,
                     tooltip = tooltip, dynamicTicks = dynamicTicks,
                     layerData = layerData, originalData = originalData,
                     source = source)
}

#' @export
as.data.frame.ArrayTimeBySpecies <- function(x, row.names = NULL,
                                             optional = FALSE, ...) {
    t <- as.numeric(rownames(x))
    if (any(is.na(t))) {
        t <- seq_len(nrow(x))
    }
    sp_names <- colnames(x)
    mat <- unclass(x)
    data.frame(
        time = rep(t, times = ncol(mat)),
        value = c(mat),
        Species = rep(sp_names, each = nrow(mat)),
        stringsAsFactors = FALSE
    )
}

#' @export
`[.ArrayTimeBySpecies` <- function(x, i, j, ..., drop = TRUE) {
    result <- NextMethod()
    # Preserve class only if result is still a 2D matrix
    if (is.matrix(result) && length(dim(result)) == 2) {
        attr(result, "value_name") <- attr(x, "value_name")
        attr(result, "units") <- attr(x, "units")
        attr(result, "params") <- attr(x, "params")
        class(result) <- c("ArrayTimeBySpecies", "matrix", "array")
    }
    result
}

#' @export
Ops.ArrayTimeBySpecies <- function(e1, e2) {
    # Strip ArrayTimeBySpecies class so that arithmetic returns a plain matrix.
    if (is.ArrayTimeBySpecies(e1)) e1 <- unclass_time(e1)
    if (!missing(e2) && is.ArrayTimeBySpecies(e2)) e2 <- unclass_time(e2)
    op <- match.fun(.Generic)
    if (missing(e2)) op(e1) else op(e1, e2)
}

# Helper to strip all ArrayTimeBySpecies attributes
unclass_time <- function(x) {
    x <- unclass(x)
    attr(x, "value_name") <- NULL
    attr(x, "units") <- NULL
    attr(x, "params") <- NULL
    x
}
