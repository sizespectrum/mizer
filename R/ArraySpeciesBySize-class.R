# ArraySpeciesBySize S3 class for species x size arrays
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

#' S3 class for species x size rate arrays
#'
#' Many functions in mizer return two-dimensional arrays (species x size)
#' holding rates like encounter rate, feeding level, growth rate, mortality etc.
#' The `ArraySpeciesBySize` class wraps these arrays to provide convenient
#' `print()`, `summary()`, `plot()`, and `as.data.frame()` methods.
#'
#' A `ArraySpeciesBySize` object behaves just like a regular matrix for
#' arithmetic operations and subsetting. It carries two lightweight attributes:
#' \itemize{
#'   \item `value_name` – a human-readable name for the value
#'       (e.g. "Encounter rate").
#'   \item `units` – the units of the rate (e.g. "g/year").
#' }
#'
#' @param x A matrix (species x size).
#' @param value_name A string giving the human-readable name for the value.
#' @param units A string giving the units (e.g. "g/year", "1/year").
#' @param params A `MizerParams` object. Currently unused but reserved for
#'   future extensions.
#'
#' @return A `ArraySpeciesBySize` object (inherits from `matrix` and `array`).
#' @export
#' @examples
#' \donttest{
#' enc <- getEncounter(NS_params)
#' is.ArraySpeciesBySize(enc)
#' summary(enc)
#' }
ArraySpeciesBySize <- function(x, value_name = NULL, units = NULL,
                               params = NULL) {
    if (!is.matrix(x)) {
        stop("`x` must be a matrix.")
    }
    if (!is.null(params)) {
        dimnames(x) <- dimnames(params@metab)
    }
    structure(x,
        class = c("ArraySpeciesBySize", "matrix", "array"),
        value_name = value_name,
        units = units
    )
}

#' Test if an object is a ArraySpeciesBySize
#'
#' @param x An object to test.
#' @return `TRUE` if `x` is a `ArraySpeciesBySize` object, `FALSE` otherwise.
#' @export
#' @examples
#' is.ArraySpeciesBySize(getEncounter(NS_params))
#' is.ArraySpeciesBySize(matrix(1:4, nrow = 2))
is.ArraySpeciesBySize <- function(x) {
    inherits(x, "ArraySpeciesBySize")
}

#' @export
print.ArraySpeciesBySize <- function(x, ...) {
    value_name <- attr(x, "value_name") %||% "ArraySpeciesBySize"
    units_str <- attr(x, "units")
    dims <- dim(x)
    header <- paste0(value_name, " (", dims[1], " species x ", dims[2],
                     " sizes)")
    if (!is.null(units_str)) {
        header <- paste0(header, " [", units_str, "]")
    }
    cat(header, "\n")
    # Print a compact summary per species
    sp_names <- rownames(x)
    if (!is.null(sp_names)) {
        vals <- apply(unclass(x), 1, function(row) {
            row <- row[is.finite(row)]
            if (length(row) == 0) return("all NA/Inf")
            paste0("min=", signif(min(row), 3),
                   " mean=", signif(mean(row), 3),
                   " max=", signif(max(row), 3))
        })
        for (i in seq_along(sp_names)) {
            cat("  ", sp_names[i], ": ", vals[i], "\n", sep = "")
        }
    }
    invisible(x)
}

#' @export
summary.ArraySpeciesBySize <- function(object, ...) {
    value_name <- attr(object, "value_name") %||% "ArraySpeciesBySize"
    units_str <- attr(object, "units")
    sp_names <- rownames(object)
    mat <- unclass(object)
    
    df <- data.frame(
        Species = sp_names,
        Min = apply(mat, 1, min, na.rm = TRUE),
        Mean = apply(mat, 1, mean, na.rm = TRUE),
        Max = apply(mat, 1, max, na.rm = TRUE),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    
    result <- list(
        value_name = value_name,
        units = units_str,
        dims = dim(object),
        per_species = df
    )
    class(result) <- "summary.ArraySpeciesBySize"
    result
}

#' @export
print.summary.ArraySpeciesBySize <- function(x, ...) {
    header <- x$value_name
    if (!is.null(x$units)) {
        header <- paste0(header, " [", x$units, "]")
    }
    cat(header, "\n")
    cat(x$dims[1], "species x", x$dims[2], "sizes\n\n")
    print(x$per_species, row.names = FALSE)
    invisible(x)
}

#' Plot a ArraySpeciesBySize object
#'
#' Plots the value against size for each species, using species colours and
#' linetypes from the MizerParams object.
#'
#' @param x A `ArraySpeciesBySize` object.
#' @param params A `MizerParams` object. Used for species colours, linetypes,
#'   and size ranges. If `NULL`, a basic plot is produced.
#' @param species Character vector of species to include. `NULL` (default) means
#'   all species.
#' @param all.sizes If `FALSE` (default), values outside a species' size range
#'   (`w_min` to `w_max`) are removed. Only effective when `params` is provided.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param return_data If `TRUE`, return the data frame instead of the plot.
#' @param log_x If `TRUE` (default), use a log10 x-axis.
#' @param ... Further arguments (currently unused).
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the variables 'w', 'value', 'Species' is returned.
#' @export
#' @examples
#' \donttest{
#' plot(getEncounter(NS_params), NS_params)
#' plot(getFeedingLevel(NS_params), NS_params,
#'      species = c("Cod", "Herring"))
#' }
plot.ArraySpeciesBySize <- function(x, params = NULL, species = NULL, 
                            all.sizes = FALSE, highlight = NULL, 
                            return_data = FALSE, log_x = TRUE, ...) {
    value_name <- attr(x, "value_name") %||% "Rate"
    units_str <- attr(x, "units")
    
    # Get w grid from params or from dimnames
    if (!is.null(params)) {
        w <- params@w
        linecolour <- params@linecolour
        linetype <- params@linetype
        species_params <- params@species_params
    } else {
        w <- as.numeric(colnames(x))
        if (any(is.na(w))) w <- seq_len(ncol(x))
        linecolour <- NULL
        linetype <- NULL
        species_params <- NULL
    }
    
    all_species <- rownames(x)
    if (is.null(species)) {
        species <- all_species
    } else {
        species <- intersect(species, all_species)
        if (length(species) == 0) {
            stop("None of the selected species are in the rate array.")
        }
    }
    
    sel <- all_species %in% species
    mat <- unclass(x)[sel, , drop = FALSE]
    
    plot_dat <- data.frame(
        w = rep(w, each = sum(sel)),
        value = c(mat),
        Species = rownames(mat)
    )
    
    if (!all.sizes && !is.null(species_params)) {
        for (sp in species) {
            if (sp %in% species_params$species) {
                sp_row <- species_params[species_params$species == sp, ]
                w_min_sp <- sp_row$w_min[1]
                w_max_sp <- sp_row$w_max[1]
                plot_dat$value[plot_dat$Species == sp &
                                   (plot_dat$w < w_min_sp |
                                        plot_dat$w > w_max_sp)] <- NA
            }
        }
        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }
    
    if (return_data) return(plot_dat)
    
    y_label <- value_name
    if (!is.null(units_str)) {
        y_label <- paste0(value_name, " [", units_str, "]")
    }
    
    # Set up species as ordered factor for legend
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
            x = .data[["w"]],
            y = .data[["value"]],
            colour = .data[["Species"]],
            linetype = .data[["Species"]],
            linewidth = .data[["Species"]]
        )) +
        ggplot2::scale_y_continuous(name = y_label)
    
    if (log_x) {
        p <- p + ggplot2::scale_x_continuous(name = "Size [g]",
                                              trans = "log10")
    } else {
        p <- p + ggplot2::scale_x_continuous(name = "Size [g]")
    }
    
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
as.data.frame.ArraySpeciesBySize <- function(x, row.names = NULL,
                                     optional = FALSE, ...) {
    w <- as.numeric(colnames(x))
    if (any(is.na(w))) {
        w <- seq_len(ncol(x))
    }
    sp_names <- rownames(x)
    mat <- unclass(x)
    data.frame(
        w = rep(w, each = nrow(mat)),
        value = c(mat),
        Species = sp_names,
        stringsAsFactors = FALSE
    )
}

#' @export
`[.ArraySpeciesBySize` <- function(x, i, j, ..., drop = TRUE) {
    result <- NextMethod()
    # Preserve class only if result is still a 2D matrix
    if (is.matrix(result) && length(dim(result)) == 2) {
        attr(result, "value_name") <- attr(x, "value_name")
        attr(result, "units") <- attr(x, "units")
        class(result) <- c("ArraySpeciesBySize", "matrix", "array")
    }
    result
}

#' @export
Ops.ArraySpeciesBySize <- function(e1, e2) {
    # Strip ArraySpeciesBySize class so that arithmetic returns a plain matrix.
    # We unclass both operands and call the generic directly.
    if (is.ArraySpeciesBySize(e1)) e1 <- unclass_rate(e1)
    if (!missing(e2) && is.ArraySpeciesBySize(e2)) e2 <- unclass_rate(e2)
    op <- match.fun(.Generic)
    if (missing(e2)) op(e1) else op(e1, e2)
}

# Helper to strip all ArraySpeciesBySize attributes
unclass_rate <- function(x) {
    x <- unclass(x)
    attr(x, "value_name") <- NULL
    attr(x, "units") <- NULL
    x
}
