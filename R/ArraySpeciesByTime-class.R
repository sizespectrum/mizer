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
#' A `ArraySpeciesByTime` object behaves just like a regular matrix for
#' arithmetic operations and subsetting. It carries two lightweight attributes:
#' \itemize{
#'   \item `value_name` – a human-readable name for the value
#'       (e.g. "Biomass").
#'   \item `units` – the units of the value (e.g. "g", "g/year").
#' }
#'
#' @param x A matrix (time x species).
#' @param value_name A string giving the human-readable name for the value.
#' @param units A string giving the units (e.g. "g", "g/year").
#' @param sim A `MizerSim` object. Currently unused but reserved for
#'   future extensions.
#'
#' @return A `ArraySpeciesByTime` object (inherits from `matrix` and `array`).
#' @export
#' @examples
#' \donttest{
#' bio <- getBiomass(NS_sim)
#' is.ArraySpeciesByTime(bio)
#' summary(bio)
#' }
ArraySpeciesByTime <- function(x, value_name = NULL, units = NULL,
                               sim = NULL) {
    if (!is.matrix(x)) {
        stop("`x` must be a matrix.")
    }
    structure(x,
        class = c("ArraySpeciesByTime", "matrix", "array"),
        value_name = value_name,
        units = units
    )
}

#' Test if an object is a ArraySpeciesByTime
#'
#' @param x An object to test.
#' @return `TRUE` if `x` is a `ArraySpeciesByTime` object, `FALSE` otherwise.
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

#' Plot a ArraySpeciesByTime object
#'
#' Plots the value against time for each species, using species colours and
#' linetypes from the MizerParams object.
#'
#' @param x A `ArraySpeciesByTime` object.
#' @param params A `MizerParams` object. Used for species colours and
#'   linetypes. If `NULL`, a basic plot is produced.
#' @param species Character vector of species to include. `NULL` (default) means
#'   all species.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param return_data If `TRUE`, return the data frame instead of the plot.
#' @param ... Further arguments (currently unused).
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with columns 'time', 'value', 'Species' is returned.
#' @export
#' @examples
#' \donttest{
#' plot(getBiomass(NS_sim), NS_params)
#' plot(getYield(NS_sim), NS_params, species = c("Cod", "Herring"))
#' }
plot.ArraySpeciesByTime <- function(x, params = NULL, species = NULL,
                                    highlight = NULL, return_data = FALSE,
                                    ...) {
    value_name <- attr(x, "value_name") %||% "Value"
    units_str <- attr(x, "units")

    t <- as.numeric(rownames(x))
    if (any(is.na(t))) t <- seq_len(nrow(x))

    if (!is.null(params)) {
        linecolour <- params@linecolour
        linetype <- params@linetype
    } else {
        linecolour <- NULL
        linetype <- NULL
    }

    all_species <- colnames(x)
    if (is.null(species)) {
        species <- all_species
    } else {
        species <- intersect(species, all_species)
        if (length(species) == 0) {
            stop("None of the selected species are in the array.")
        }
    }

    sel <- all_species %in% species
    mat <- unclass(x)[, sel, drop = FALSE]

    plot_dat <- data.frame(
        time = rep(t, times = sum(sel)),
        value = c(mat),
        Species = rep(colnames(mat), each = nrow(mat)),
        stringsAsFactors = FALSE
    )

    if (return_data) return(plot_dat)

    y_label <- value_name
    if (!is.null(units_str) && nzchar(units_str)) {
        y_label <- paste0(value_name, " [", units_str, "]")
    }

    if (!is.null(linecolour)) {
        legend_levels <- intersect(names(linecolour), plot_dat$Species)
    } else {
        legend_levels <- unique(plot_dat$Species)
    }
    plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)

    linesize <- rep(0.8, length(legend_levels))
    names(linesize) <- legend_levels
    if (!is.null(highlight)) {
        linesize[highlight] <- 1.6
    }

    p <- ggplot2::ggplot(plot_dat, ggplot2::aes(group = .data[["Species"]])) +
        ggplot2::geom_line(ggplot2::aes(
            x = .data[["time"]],
            y = .data[["value"]],
            colour = .data[["Species"]],
            linetype = .data[["Species"]],
            linewidth = .data[["Species"]]
        )) +
        ggplot2::scale_x_continuous(name = "Time [years]") +
        ggplot2::scale_y_continuous(name = y_label)

    if (!is.null(linecolour)) {
        p <- p +
            ggplot2::scale_colour_manual(values = linecolour[legend_levels])
    }
    if (!is.null(linetype)) {
        p <- p +
            ggplot2::scale_linetype_manual(values = linetype[legend_levels])
    }
    p <- p + ggplot2::scale_discrete_manual("linewidth", values = linesize)

    p
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
    x
}
