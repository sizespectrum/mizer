#' Extract the `major.minor` part of a version string
#'
#' Returns `""` for anything that does not start with `major.minor`, so that a
#' missing or malformed record simply counts as "not the current series".
#'
#' @param version A character vector of version strings.
#' @noRd
minor_version <- function(version) {
    ifelse(grepl("^[0-9]+\\.[0-9]+", version),
           sub("^([0-9]+\\.[0-9]+).*$", "\\1", version),
           "")
}

#' Show a startup message the first time a new mizer series is loaded
#'
#' Compares the installed version against the version last recorded in the
#' user's mizer config directory (`tools::R_user_dir("mizer", "config")`) and
#' prints a startup message and updates the record if they differ. Only the
#' `major.minor` part is compared, because the message links to the release
#' announcement for the series and there is one of those per minor release.
#' A patch release or a development version therefore does not re-show an
#' announcement the user has already seen. Any failure (e.g. an unwritable
#' config directory) is swallowed so that `library(mizer)` can never be broken
#' by this.
#'
#' @param libname unused
#' @param pkgname unused
#' @noRd
.onAttach <- function(libname, pkgname) {
    tryCatch({
        current <- as.character(utils::packageVersion("mizer"))
        config_dir <- tools::R_user_dir("mizer", "config")
        version_file <- file.path(config_dir, "last_seen_version")
        last_seen <- if (file.exists(version_file)) {
            readLines(version_file, warn = FALSE)
        } else {
            ""
        }
        if (!identical(minor_version(last_seen), minor_version(current))) {
            packageStartupMessage(
                "This is mizer version ", current,
                ". See https://blog.mizer.sizespectrum.org/posts/2026-08-30-mizer-3-4-announcement/ for what has changed."
            )
            dir.create(config_dir, showWarnings = FALSE, recursive = TRUE)
            # Use write() rather than writeLines()/cat() so that R CMD check's
            # startup-function heuristic does not mistake this file write for a
            # console message (it flags those calls by name, not by argument).
            write(current, file = version_file)
        }
    },
    warning = function(w) invisible(NULL),
    error = function(e) invisible(NULL))
}
