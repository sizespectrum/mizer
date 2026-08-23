# Helpers for the articles whose results are pre-computed
# (dev_scripts/save_dynamic_stability.R and
# dev_scripts/save_discontinuous_rates.R).
#
# Those scripts hold the same code as their article, which shows it with
# `eval = FALSE` and displays the saved objects instead of running them. The
# saved file is therefore the only place where a change in mizer's behaviour
# becomes visible, and it is easy to overwrite it without noticing that a number
# the article quotes in its text has moved. `save_with_report()` compares the
# newly computed results against the ones already on disk and says what changed
# before writing.

# Slots that change on every rebuild and say nothing about behaviour.
.ignored_slots <- c("time_created", "time_modified")

# Largest relative difference between two numeric vectors, treating 0 vs 0 as 0.
# NA, NaN and infinite entries are compared exactly: a value that was Inf and is
# now finite, or an NA that has moved, counts as a difference of its own rather
# than producing a NaN that would silently escape the comparison.
.rel_diff <- function(old, new) {
    special_old <- !is.finite(old)
    special_new <- !is.finite(new)
    if (!identical(special_old, special_new)) return(Inf)
    if (any(special_old) &&
            !identical(old[special_old], new[special_old])) {
        return(Inf)
    }
    keep <- !special_old
    if (!any(keep)) return(0)
    old <- old[keep]
    new <- new[keep]
    scale <- pmax(abs(old), abs(new))
    d <- ifelse(scale == 0, 0, abs(new - old) / scale)
    max(d)
}

.fmt <- function(x) format(signif(x, 4), trim = TRUE)

.describe_dim <- function(x) {
    if (!is.null(dim(x))) paste(dim(x), collapse = " x ") else paste(length(x))
}

#' Compare two saved result objects
#'
#' @param old,new The objects to compare.
#' @param tol Relative differences at or below this are not reported.
#' @param path The name of the object, used to build the reported labels.
#' @param depth Recursion depth, used to stop on an unexpectedly nested object.
#' @return A character vector with one entry per difference found, empty if
#'   the two agree to within `tol`.
#' @noRd
compare_saved <- function(old, new, tol = 1e-8, path = "", depth = 0L) {
    label <- if (nzchar(path)) path else "result"

    # A structure deeper than anything these articles save. Rather than risk
    # running out of stack, fall back to a yes/no answer.
    if (depth > 8L) {
        if (!isTRUE(all.equal(old, new, tolerance = tol))) {
            return(paste0(label, ": differs (not examined further)"))
        }
        return(character(0))
    }

    if (is.null(old) != is.null(new)) {
        return(paste0(label, ": ", if (is.null(old)) "added" else "removed"))
    }
    if (is.null(old)) return(character(0))

    if (!identical(class(old), class(new))) {
        return(paste0(label, ": class changed, ",
                      paste(class(old), collapse = "/"), " -> ",
                      paste(class(new), collapse = "/")))
    }

    # A MizerParams object holds far more than an article ever shows, and
    # walking all of it would bury the report. What matters for an article is
    # the state the model is in and the parameters behind it.
    if (is(old, "MizerParams")) {
        parts <- list(initialN = initialN,
                      initialNResource = initialNResource,
                      initial_effort = initial_effort,
                      species_params = species_params)
        return(unlist(lapply(names(parts), function(nm) {
            compare_saved(parts[[nm]](old), parts[[nm]](new), tol = tol,
                          path = paste0(label, "@", nm), depth = depth + 1L)
        }), use.names = FALSE))
    }

    if (isS4(old)) {
        slots <- setdiff(methods::slotNames(old), .ignored_slots)
        return(unlist(lapply(slots, function(s) {
            compare_saved(methods::slot(old, s), methods::slot(new, s),
                          tol = tol, path = paste0(label, "@", s),
                          depth = depth + 1L)
        }), use.names = FALSE))
    }

    if (is.function(old)) {
        if (!identical(deparse(old), deparse(new))) {
            return(paste0(label, ": function definition changed"))
        }
        return(character(0))
    }

    if (is.list(old)) {
        diffs <- character(0)
        if (is.data.frame(old) && nrow(old) != nrow(new)) {
            diffs <- c(diffs, paste0(label, ": ", nrow(old), " rows -> ",
                                     nrow(new), " rows"))
        }
        nms_old <- names(old)
        nms_new <- names(new)
        if (is.null(nms_old) || is.null(nms_new)) {
            if (length(old) != length(new)) {
                return(paste0(label, ": length ", length(old), " -> ",
                              length(new)))
            }
            nms_old <- nms_new <- seq_along(old)
            get <- function(x, i) x[[i]]
        } else {
            gone <- setdiff(nms_old, nms_new)
            added <- setdiff(nms_new, nms_old)
            if (length(gone)) {
                diffs <- c(diffs, paste0(label, ": removed ",
                                         paste(gone, collapse = ", ")))
            }
            if (length(added)) {
                diffs <- c(diffs, paste0(label, ": added ",
                                         paste(added, collapse = ", ")))
            }
            get <- function(x, i) x[[i]]
        }
        common <- intersect(nms_old, nms_new)
        sep <- if (nzchar(path)) "$" else ""
        for (nm in common) {
            diffs <- c(diffs,
                       compare_saved(get(old, nm), get(new, nm), tol = tol,
                                     path = paste0(path, sep, nm),
                                     depth = depth + 1L))
        }
        return(diffs)
    }

    # From here on we compare values only. Attributes travel with mizer's
    # arrays (a `params` object, units, a `comment` marking a frozen rate) and
    # are not what a reader of the article sees.
    if (!identical(dim(old), dim(new)) || length(old) != length(new)) {
        return(paste0(label, ": size ", .describe_dim(old), " -> ",
                      .describe_dim(new)))
    }

    if (is.numeric(old) || is.complex(old)) {
        d <- .rel_diff(as.vector(old), as.vector(new))
        if (d <= tol) return(character(0))
        if (length(old) == 1) {
            return(paste0(label, ": ", .fmt(old), " -> ", .fmt(new),
                          "  (rel. diff ", .fmt(d), ")"))
        }
        n <- sum(as.vector(old) != as.vector(new), na.rm = TRUE)
        return(paste0(label, ": ", n, " of ", length(old),
                      " values differ, largest rel. diff ", .fmt(d)))
    }

    # Character, logical, factor and anything else.
    if (length(old) == 0) return(character(0))
    same <- mapply(identical, as.vector(old), as.vector(new))
    if (all(same)) return(character(0))
    which_diff <- which(!same)
    example <- which_diff[[1]]
    detail <- paste0(" (e.g. [", example, "] ",
                     format(as.vector(old)[[example]]), " -> ",
                     format(as.vector(new)[[example]]), ")")
    if (length(old) == 1) {
        return(paste0(label, ": ", format(old), " -> ", format(new)))
    }
    paste0(label, ": ", length(which_diff), " of ", length(old),
           " entries differ", detail)
}

#' Save pre-computed article results, reporting what changed
#'
#' Compares the new results against the file already at `path`, reports the
#' differences, and then writes the new ones.
#'
#' @param results A named list of the objects the article displays.
#' @param path The file to write.
#' @param tol Relative differences at or below this are not reported.
#' @param max_lines The most difference lines to print.
#' @noRd
save_with_report <- function(results, path, tol = 1e-8, max_lines = 40) {
    if (file.exists(path)) {
        diffs <- compare_saved(readRDS(path), results, tol = tol)
        if (length(diffs) == 0) {
            message("No change: the recalculated results agree with ", path,
                    " to within a relative ", .fmt(tol), ".")
        } else {
            message("Changes from the previous ", path, ":")
            shown <- utils::head(diffs, max_lines)
            message(paste0("  ", shown, collapse = "\n"))
            if (length(diffs) > max_lines) {
                message("  ... and ", length(diffs) - max_lines, " more.")
            }
            message("The article quotes some of these numbers in its text. ",
                    "Re-knit it and check the inline values.")
        }
    } else {
        message("No existing ", path, ", so there is nothing to compare with.")
    }
    saveRDS(results, path)
    message("Saved ", path)
}
