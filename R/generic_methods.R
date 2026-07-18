# Documentation for methods of base generics
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

# print ----
#' Print mizer objects
#'
#' Mizer supplies `print()` methods for the array-like objects returned by many
#' rate and summary functions. These methods print a compact preview of the
#' underlying matrix, array or vector: a header reporting the value name,
#' dimensions and units, followed by the actual values, truncated to fit the
#' console when the array is large. Species are truncated to a leading
#' subset, sizes to an evenly log-spaced sample spanning the full size range
#' (because size grids are uniform in log-space, this shows small, medium and
#' large individuals rather than just the smallest), and time series to a
#' representative sample of time steps that always includes the first and
#' last. A trailing note reports how much was omitted. A three-dimensional
#' [ArrayTimeBySpeciesBySize()] object is previewed via its final time slice,
#' matching the default behaviour of [plot()] for that class.
#'
#' For full numeric access, use the object itself as an ordinary matrix, array
#' or vector, or convert it to a long data frame with [as.data.frame()].
#'
#' @param x The object to print.
#' @param ... Further arguments. They are currently ignored by the mizer
#'   methods.
#'
#' @return The printed object, invisibly.
#'
#' @usage print(x, ...)
#' @aliases print.ArraySpeciesBySize print.ArrayTimeBySpecies print.ArrayTimeBySpeciesBySize print.summary.ArraySpeciesBySize print.summary.ArrayTimeBySpecies print.summary.ArrayTimeBySpeciesBySize
#' @seealso [summary()], [as.data.frame()], [plot()],
#'   [ArraySpeciesBySize()], [ArrayTimeBySpecies()],
#'   [ArrayTimeBySpeciesBySize()]
#'
#' @examples
#' \donttest{
#' enc <- getEncounter(NS_params)
#' print(enc)
#'
#' biomass <- getBiomass(NS_sim)
#' print(biomass)
#' }
#' @name print
NULL

# summary ----
#' Summarise mizer objects
#'
#' Mizer provides `summary()` methods for model objects and for the specialised
#' array classes returned by many mizer functions.
#'
#' For a [MizerParams()] object, `summary()` prints the model metadata, size
#' grids, selected species parameters and fishing gear details. For a
#' [MizerSim()] object, it first prints the parameter summary and then reports
#' the simulated time period and output interval.
#'
#' For [ArraySpeciesBySize()], [ArrayTimeBySpecies()] and
#' [ArrayTimeBySpeciesBySize()] objects, `summary()` returns a small list with
#' the value name, units, dimensions and a per-species data frame containing
#' minimum, mean and maximum values. Printing that summary object gives the same
#' compact table in a human-readable form.
#'
#' @param object The object to summarise.
#' @param ... Further arguments. They are currently ignored by the mizer
#'   methods.
#'
#' @return
#' For [MizerParams()] and [MizerSim()], the object is returned invisibly.
#' For array objects, a list of class `summary.ArraySpeciesBySize`,
#' `summary.ArrayTimeBySpecies` or `summary.ArrayTimeBySpeciesBySize`.
#'
#' @usage summary(object, ...)
#' @aliases summary.ArraySpeciesBySize summary.ArrayTimeBySpecies summary.ArrayTimeBySpeciesBySize summary.MizerSim summary.MizerParams
#' @seealso [print()], [as.data.frame()], [str()], [MizerParams()], [MizerSim()],
#'   [ArraySpeciesBySize()], [ArrayTimeBySpecies()],
#'   [ArrayTimeBySpeciesBySize()]
#'
#' @examples
#' \donttest{
#' summary(NS_params)
#' summary(NS_sim)
#' summary(getEncounter(NS_params))
#' summary(getFMort(NS_sim))
#' }
#' @name summary
NULL

# as.data.frame ----
#' Convert mizer arrays to data frames
#'
#' The `as.data.frame()` methods for mizer array classes turn matrix- and
#' array-like results into tidy long-form data frames, with one row per
#' observed combination of species, size and/or time. The numeric result is
#' always stored in a column called `value`.
#'
#' The returned columns are:
#' \itemize{
#'   \item `ArraySpeciesBySize`: `w`, `value`, `Species`.
#'   \item `ArrayTimeBySpecies`: `time`, `value`, `Species`.
#'   \item `ArrayTimeBySpeciesBySize`: `time`, `Species`, `w`, `value`.
#' }
#'
#' If the original object has non-numeric or missing dimension names, sequential
#' indices are used for the `time` or `w` columns. Species names are taken from
#' the row, column or dimension names of the original object.
#'
#' @param x An `ArraySpeciesBySize`, `ArrayTimeBySpecies` or
#'   `ArrayTimeBySpeciesBySize` object.
#' @param row.names Optional and included only for compatibility with the base
#'   generic. `NULL` or a character vector giving the row names for the
#'   data frame.
#' @param optional Optional and included only for compatibility with the base
#'   generic. A logical value. If `TRUE`, setting row names and converting
#'   column names (to syntactic names) is optional.
#' @param ... Further arguments. They are currently ignored by the mizer
#'   methods.
#'
#' @return A data frame in long format.
#'
#' @usage as.data.frame(x, row.names = NULL, optional = FALSE, ...)
#' @aliases as.data.frame.ArraySpeciesBySize as.data.frame.ArrayTimeBySpecies as.data.frame.ArrayTimeBySpeciesBySize
#' @seealso [print()], [summary()], [str()], [plot()], [ArraySpeciesBySize()],
#'   [ArrayTimeBySpecies()], [ArrayTimeBySpeciesBySize()]
#'
#' @examples
#' \donttest{
#' enc <- getEncounter(NS_params)
#' head(as.data.frame(enc))
#'
#' biomass <- getBiomass(NS_sim)
#' head(as.data.frame(biomass))
#' }
#' @name as.data.frame
NULL

# str ----
#' Display the structure of mizer objects
#'
#' Mizer provides `str()` methods for [MizerParams()] and [MizerSim()] objects,
#' as well as [ArraySpeciesBySize()], [ArrayTimeBySpecies()] and [ArrayTimeBySpeciesBySize()]
#' objects. These methods produce a clean, compact overview of the object's structure
#' without polluting the console with large amounts of internal data.
#'
#' @param object The object to display the structure of.
#' @param max.level Maximum level of nesting to print. Defaults to `NA` (no limit).
#' @param ... Further arguments. They are passed to the default `str()` method.
#'
#' @return `NULL`, invisibly.
#'
#' @usage str(object, max.level = NA, ...)
#' @aliases str.ArraySpeciesBySize str.ArrayTimeBySpecies str.ArrayTimeBySpeciesBySize str.MizerSim str.MizerParams
#' @seealso [print()], [as.data.frame()], [summary()], [plot()], [MizerParams()], [MizerSim()], [ArraySpeciesBySize()],
#'   [ArrayTimeBySpecies()], [ArrayTimeBySpeciesBySize()]
#'
#' @examples
#' \donttest{
#' str(NS_params)
#' str(NS_sim)
#' str(getEncounter(NS_params))
#' }
#' @name str
NULL

# Print truncation helpers ---------------------------------------------------
#
# Shared, non-exported machinery used by the print.Array*() methods (defined
# in the various *-class.R files) to keep their console output compact for
# large models or long simulations.
mizer_print_defaults <- list(
    time_max = 8L,
    time_threshold = 12L,
    species_head = 8L,
    species_threshold = 13L,
    size_max = 10L,
    size_min = 2L
)

# Choose evenly spaced indices out of 1:n, always including 1 and n, unless
# n is at most `threshold` (which defaults to max_n itself). Because mizer's
# size grids are uniform in log-space, spacing the *indices* evenly
# reproduces an even spread of *sizes* from smallest to largest; for a
# regularly-saved time grid it reproduces an even spread of *times*.
# `threshold` lets a caller trigger truncation later than `max_n` alone would
# (e.g. only truncate once well past the display count), mirroring
# pick_head_indices()'s `threshold` argument.
pick_log_spaced_indices <- function(n, max_n, threshold = max_n) {
    if (n <= threshold) {
        return(seq_len(n))
    }
    sort(unique(round(seq(1, n, length.out = max_n))))
}

# Choose a plain leading subset of 1:n, for dimensions (like species) that
# are rarely large enough to need truncation at all.
pick_head_indices <- function(n, head_n, threshold = head_n) {
    if (n <= threshold) {
        return(seq_len(n))
    }
    seq_len(head_n)
}

# Estimate the number of characters needed to print a numeric matrix (with
# its row names) in one unwrapped block, using the same format() logic that
# base print.default() uses internally to lay out columns.
matrix_display_width <- function(mat) {
    formatted <- format(mat)
    col_widths <- apply(formatted, 2, function(col) max(nchar(col)))
    row_label_width <- max(nchar(rownames(mat)), 0)
    row_label_width + sum(col_widths + 1)
}

# Estimate the number of characters needed to print a named numeric vector
# in one unwrapped block.
vector_display_width <- function(vec) {
    formatted <- format(vec)
    elem_widths <- pmax(nchar(names(vec) %||% ""), nchar(formatted))
    sum(elem_widths + 1)
}

# Choose how many log-spaced columns/elements (out of n) fit within the
# console width, by trying counts from max_k down to min_k and asking
# `width_fn(k)` (supplied by the caller, using matrix_display_width() or
# vector_display_width() on the candidate subset) how wide that would print.
fit_log_spaced_k <- function(n, max_k, min_k, width_fn,
                             width = getOption("width", 80L)) {
    k_start <- min(max_k, n)
    if (k_start <= min_k) {
        return(k_start)
    }
    for (k in k_start:min_k) {
        if (width_fn(k) <= width) {
            return(k)
        }
    }
    min_k
}

# Parse dimnames labels as numeric, falling back to a sequential index when
# they are missing or not parseable (mirrors the pattern already used by e.g.
# get_ArrayResourceBySize_w() and get_ArraySpeciesBySize_w()).
parse_numeric_labels <- function(labels, n) {
    w <- suppressWarnings(as.numeric(labels))
    if (is.null(labels) || anyNA(w)) {
        w <- seq_len(n)
    }
    w
}

# Format the trailing note reporting how much of an array was omitted from a
# print preview, e.g. "... showing 10 of 226 sizes (0.001-4e+04 g,
# log-spaced); use as.data.frame() for the full data."
format_truncation_note <- function(shown, total, label, detail = NULL) {
    note <- paste0("... showing ", shown, " of ", total, " ", label)
    if (!is.null(detail) && nzchar(detail)) {
        note <- paste0(note, " (", detail, ")")
    }
    paste0(note, "; use as.data.frame() for the full data.")
}

# Format the size-range detail used inside format_truncation_note() for the
# size dimension. Weights are always in grams (see mizer's units convention).
format_size_range_detail <- function(w) {
    paste0(signif(min(w), 3), "-", signif(max(w), 3), " g, log-spaced")
}

# Format the time-range detail used inside format_truncation_note() for an
# evenly spaced sample of time steps.
format_time_range_detail <- function(times) {
    paste0(format(min(times), trim = TRUE), "-", format(max(times), trim = TRUE),
          ", evenly spaced")
}

