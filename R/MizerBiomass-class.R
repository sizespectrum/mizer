# Metadata carried by explicitly selected biomass measurements. The ordinary
# species-only getter retains its existing numeric-vector and array interface.

measurement_attributes <- function() {
    c("series", "coverage", "unavailable")
}

subset_measurement_metadata <- function(result, original, ids, times = NULL) {
    series <- attr(original, "series")
    if (is.null(series)) return(result)
    series <- series[match(ids, series$id), , drop = FALSE]
    rownames(series) <- NULL
    coverage <- attr(original, "coverage")
    coverage <- coverage[coverage$id %in% ids, , drop = FALSE]
    if (!is.null(times)) {
        coverage <- coverage[as.character(coverage$time) %in%
                                 as.character(times), , drop = FALSE]
    }
    attr(result, "series") <- series
    attr(result, "coverage") <- coverage
    attr(result, "unavailable") <- coverage[!coverage$available, , drop = FALSE]
    result
}

measurement_caption <- function(x) {
    unavailable <- attr(x, "unavailable")
    if (is.null(unavailable) || !nrow(unavailable)) return(NULL)
    paste("Unavailable:", paste(unique(paste0(
        unavailable$id, " (", unavailable$reason, ")"
    )), collapse = "; "))
}

check_measurement_total <- function(x) {
    coverage <- attr(x, "coverage")
    if (!is.null(coverage) && (any(!coverage$available) || anyNA(x))) {
        stop("A total is unavailable because requested biomass series are ",
             "incomplete. Inspect the measurement's coverage metadata.",
             call. = FALSE)
    }
    ids <- attr(x, "series")$source
    if (!is.null(ids) && anyDuplicated(ids)) {
        stop("A biomass total cannot count the same physical source twice.",
             call. = FALSE)
    }
    invisible(TRUE)
}

annotate_measurement_plot <- function(plot, x) {
    caption <- measurement_caption(x)
    if (!is.null(caption)) plot <- plot + ggplot2::labs(caption = caption)
    for (attribute in measurement_attributes()) {
        attr(plot, attribute) <- attr(x, attribute)
    }
    plot
}

#' Plot a biomass measurement at one state
#'
#' Explicit resource or component inclusion in [getBiomass()] returns a named
#' numeric `MizerBiomass` vector carrying its own series labels, colours and
#' coverage. This method draws those measurements without discovering additional
#' series from the model. Unavailable measurements remain in the returned data
#' and are explained in the plot caption. A total requires complete coverage.
#'
#' @param x A `MizerBiomass` vector returned by [getBiomass()].
#' @param species Series identifiers to draw; NULL draws every measured series.
#' @param total Add the total of every measured series, before selection.
#' @param log_y Use a logarithmic biomass axis.
#' @param ylim Numeric lower and upper biomass-axis limits; NA chooses them.
#' @param return_data Return the plotting data instead of a plot.
#' @param ... Currently unused.
#' @return A mizer plot, or a data frame when `return_data = TRUE`.
#' @keywords internal
#' @export
plot.MizerBiomass <- function(x, species = NULL, total = FALSE,
                              log_y = TRUE, ylim = c(NA, NA),
                              return_data = FALSE, ...) {
    metadata <- attr(x, "series")
    ids <- names(x)
    if (!is.null(species)) {
        if (any(!species %in% ids)) stop("Unknown biomass series selection.")
        ids <- unique(species)
    }
    frame <- data.frame(Series = names(x), Biomass = as.numeric(x),
                        stringsAsFactors = FALSE)
    if (total) {
        check_measurement_total(x)
        frame <- rbind(frame, data.frame(Series = "Total", Biomass = sum(x)))
        metadata <- rbind(metadata, transform(metadata[1, , drop = FALSE],
            id = "Total", label = "Total", colour = "black", linetype = "solid"
        ))
        ids <- c(ids, "Total")
    }
    frame <- frame[frame$Series %in% ids, , drop = FALSE]
    metadata <- metadata[match(frame$Series, metadata$id), , drop = FALSE]
    frame$Label <- metadata$label
    frame$Series <- factor(frame$Series, levels = unique(frame$Series))
    for (attribute in measurement_attributes()) {
        attr(frame, attribute) <- attr(x, attribute)
    }
    if (return_data) return(frame)
    shown <- frame[is.finite(frame$Biomass) &
                       (!log_y | frame$Biomass > 0), , drop = FALSE]
    colours <- stats::setNames(metadata$colour, metadata$id)
    labels <- stats::setNames(metadata$label, metadata$id)
    p <- ggplot2::ggplot(shown, ggplot2::aes(
        x = .data[["Series"]], y = .data[["Biomass"]],
        colour = .data[["Series"]]
    )) + ggplot2::geom_point(size = 3) +
        ggplot2::scale_x_discrete(labels = labels,
                                  limits = levels(frame$Series)) +
        ggplot2::scale_colour_manual(values = colours, limits = names(colours),
                                     guide = "none") +
        ggplot2::scale_y_continuous(
            trans = if (log_y) "log10" else "identity", limits = ylim
        ) + ggplot2::labs(x = NULL, y = "Biomass [g]")
    p <- make_mizer_plot(p, c("Label", "Biomass"))
    annotate_measurement_plot(p, x)
}

#' @export
print.MizerBiomass <- function(x, ...) {
    value <- as.numeric(x)
    names(value) <- names(x)
    print(value, ...)
    caption <- measurement_caption(x)
    if (!is.null(caption)) cat(caption, "\n")
    invisible(x)
}

#' @export
`[.MizerBiomass` <- function(x, i, ...) {
    result <- NextMethod()
    for (attribute in c("params", "value_name", "units", "type")) {
        attr(result, attribute) <- attr(x, attribute)
    }
    result <- subset_measurement_metadata(result, x, names(result))
    class(result) <- c("MizerBiomass", "numeric")
    result
}

#' @export
Ops.MizerBiomass <- function(e1, e2) {
    plain <- function(x) {
        value <- as.numeric(x)
        names(value) <- names(x)
        value
    }
    if (inherits(e1, "MizerBiomass")) e1 <- plain(e1)
    if (!missing(e2) && inherits(e2, "MizerBiomass")) e2 <- plain(e2)
    op <- match.fun(.Generic) # nolint: object_usage_linter.
    if (missing(e2)) op(e1) else op(e1, e2)
}

#' @export
Math.MizerBiomass <- function(x, ...) {
    value <- as.numeric(x)
    names(value) <- names(x)
    fun <- match.fun(.Generic) # nolint: object_usage_linter.
    fun(value, ...)
}

#' @export
as.data.frame.MizerBiomass <- function(
    x, row.names = NULL, optional = FALSE, ... # nolint: object_name_linter.
) {
    result <- attr(x, "series")
    result$value <- as.numeric(x)
    coverage <- attr(x, "coverage")
    index <- match(result$id, coverage$id)
    result$available <- coverage$available[index]
    result$reason <- coverage$reason[index]
    if (!is.null(row.names)) rownames(result) <- row.names
    result
}

#' @export
str.MizerBiomass <- function(object, ...) {
    attr(object, "params") <- NULL
    utils::str(unclass(object), ...)
    invisible(NULL)
}
