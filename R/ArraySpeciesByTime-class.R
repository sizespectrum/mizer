# ArraySpeciesByTime S3 class for time x species arrays
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' S3 class for time x species arrays
#'
#' Some functions in mizer return two-dimensional arrays (time x species)
#' holding quantities like biomass, abundance, or yield rate through time.
#' The `ArraySpeciesByTime` class wraps these arrays to provide convenient
#' `print()`, `summary()`, `plot()`, and `as.data.frame()` methods.
#'
#' An `ArraySpeciesByTime` object behaves just like a regular matrix for
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
#' @return An `ArraySpeciesByTime` object (inherits from `matrix` and `array`).
#' @export
#' @examples
#' \donttest{
#' bio <- getBiomass(NS_sim)
#' is.ArraySpeciesByTime(bio)
#' summary(bio)
#' }
ArraySpeciesByTime <- function(x, value_name = NULL, units = NULL,
                               params = NULL) {
    if (!is.matrix(x)) {
        stop("`x` must be a matrix.")
    }
    structure(x,
        class = c("ArraySpeciesByTime", "matrix", "array"),
        value_name = value_name,
        units = units,
        params = params
    )
}

#' Test if an object is a ArraySpeciesByTime
#'
#' @param x An object to test.
#' @return `TRUE` if `x` is an `ArraySpeciesByTime` object, `FALSE` otherwise.
#' @export
#' @examples
#' is.ArraySpeciesByTime(getBiomass(NS_sim))
#' is.ArraySpeciesByTime(matrix(1:4, nrow = 2))
is.ArraySpeciesByTime <- function(x) {
    inherits(x, "ArraySpeciesByTime")
}

#' @export
print.ArraySpeciesByTime <- function(x, ...) {
    value_name <- attr(x, "value_name") %||% "ArraySpeciesByTime"
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
summary.ArraySpeciesByTime <- function(object, ...) {
    value_name <- attr(object, "value_name") %||% "ArraySpeciesByTime"
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
    class(result) <- "summary.ArraySpeciesByTime"
    result
}

#' @export
print.summary.ArraySpeciesByTime <- function(x, ...) {
    header <- x$value_name
    if (!is.null(x$units) && nzchar(x$units)) {
        header <- paste0(header, " [", x$units, "]")
    }
    cat(header, "\n")
    cat(x$dims[1], "times x", x$dims[2], "species\n\n")
    print(x$per_species, row.names = FALSE)
    invisible(x)
}

#' Plot an ArraySpeciesByTime object
#'
#' Plots the value against time for each species, using species colours and
#' linetypes from the MizerParams object. This method supports the same
#' arguments as [plotBiomass()], so that
#' `plot(getBiomass(sim, ...), ...)` is equivalent to
#' `plotBiomass(sim, ...)`.
#'
#' @param x An `ArraySpeciesByTime` object.
#' @param species Character vector of species to include. `NULL` (default) means
#'   all species.
#' @param start_time The first time to be plotted. Default (`NULL`) is the
#'   beginning of the time series.
#' @param end_time The last time to be plotted. Default (`NULL`) is the end of
#'   the time series.
#' @param y_ticks The approximate number of ticks desired on the y axis.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param total A boolean value that determines whether the total from all
#'   species is plotted as well. Default is FALSE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param log If `TRUE` (default), use a log10 y-axis. Appropriate for
#'   quantities like biomass, SSB, and yield that span orders of magnitude.
#'   Set to `FALSE` for a linear y-axis.
#' @param return_data If `TRUE`, return the data frame instead of the plot.
#' @param ... Further arguments (currently unused).
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with columns 'Year', value_name, 'Species', 'Legend' is returned.
#' @export
#' @examples
#' \donttest{
#' plot(getBiomass(NS_sim))
#' plot(getBiomass(NS_sim), species = c("Cod", "Herring"), total = TRUE)
#' plot(getYield(NS_sim), species = c("Cod", "Herring"))
#' }
plot.ArraySpeciesByTime <- function(x, species = NULL,
                                    start_time = NULL, end_time = NULL,
                                    y_ticks = 6, ylim = c(NA, NA),
                                    total = FALSE, background = TRUE,
                                    highlight = NULL, log = TRUE,
                                    return_data = FALSE,
                                    ...) {
    value_name <- attr(x, "value_name") %||% "Value"
    units_str <- attr(x, "units")
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
        t <- t[t <= end_time]
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

    # Construct y_label
    y_label <- value_name
    if (!is.null(units_str) && nzchar(units_str)) {
        y_label <- paste0(value_name, " [", units_str, "]")
    }

    # Implement ylim and a minimal cutoff; bring columns in desired order
    min_value <- 1e-20
    bm <- bm[bm$value >= min_value &
                 (is.na(ylim[1]) | bm$value >= ylim[1]) &
                 (is.na(ylim[2]) | bm$value <= ylim[2]), c(1, 3, 2)]
    names(bm) <- c("Year", value_name, "Species")

    # Select species
    plot_dat <- bm[bm$Species %in% c("Total", species), ]
    plot_dat$Legend <- plot_dat$Species

    if (background && !is.null(params) && anyNA(params@A)) {
        # Add background species in light grey
        bkgrd_sp <- params@species_params$species[is.na(params@A)]
        if (length(bkgrd_sp) > 0) {
            bm_bkgrd <- bm[bm$Species %in% bkgrd_sp, ]
            bm_bkgrd$Legend <- "Background"
            plot_dat <- rbind(plot_dat, bm_bkgrd)
        }
    }

    if (return_data) return(plot_dat)

    plotDataFrame(plot_dat, params, xlab = "Year", ylab = y_label,
                  ytrans = if (log) "log10" else "identity", ylim = ylim,
                  y_ticks = y_ticks, highlight = highlight,
                  legend_var = "Legend")
}

#' @rdname plotly
#' @export
#' @examples
#' \donttest{
#' plotly(getBiomass(NS_sim))
#' }
plotly.ArraySpeciesByTime <- function(x, ...) {
    plotly::ggplotly(plot(x, ...))
}

#' @export
as.data.frame.ArraySpeciesByTime <- function(x, row.names = NULL,
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
`[.ArraySpeciesByTime` <- function(x, i, j, ..., drop = TRUE) {
    result <- NextMethod()
    # Preserve class only if result is still a 2D matrix
    if (is.matrix(result) && length(dim(result)) == 2) {
        attr(result, "value_name") <- attr(x, "value_name")
        attr(result, "units") <- attr(x, "units")
        attr(result, "params") <- attr(x, "params")
        class(result) <- c("ArraySpeciesByTime", "matrix", "array")
    }
    result
}

#' @export
Ops.ArraySpeciesByTime <- function(e1, e2) {
    # Strip ArraySpeciesByTime class so that arithmetic returns a plain matrix.
    if (is.ArraySpeciesByTime(e1)) e1 <- unclass_time(e1)
    if (!missing(e2) && is.ArraySpeciesByTime(e2)) e2 <- unclass_time(e2)
    op <- match.fun(.Generic)
    if (missing(e2)) op(e1) else op(e1, e2)
}

# Helper to strip all ArraySpeciesByTime attributes
unclass_time <- function(x) {
    x <- unclass(x)
    attr(x, "value_name") <- NULL
    attr(x, "units") <- NULL
    attr(x, "params") <- NULL
    x
}
