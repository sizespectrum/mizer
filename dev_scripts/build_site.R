# Build the pkgdown website without allowing standalone vignette HTML files to
# overwrite pkgdown's formatted article pages.
#
# From the repository root, run:
#
#     Rscript dev_scripts/build_site.R
#
# Or from R:
#
#     source("dev_scripts/build_site.R")
#     build_mizer_site()

build_mizer_site <- function(pkg_root = ".") {
    pkg_root <- normalizePath(pkg_root, mustWork = TRUE)
    description <- file.path(pkg_root, "DESCRIPTION")
    vignette_dir <- file.path(pkg_root, "vignettes")

    if (!file.exists(description) || !dir.exists(vignette_dir)) {
        stop("`", pkg_root, "` is not the root of the mizer package.",
             call. = FALSE)
    }

    html_files <- list.files(vignette_dir,
                             pattern = "[.]html$",
                             all.files = TRUE,
                             full.names = TRUE,
                             recursive = FALSE)
    if (length(html_files)) {
        removed <- file.remove(html_files)
        if (any(removed)) {
            message("Removed standalone vignette HTML file(s):\n",
                    paste0("  ", basename(html_files[removed]),
                           collapse = "\n"),
                    "\nThese are generated, Git-ignored build artefacts.")
        }
        if (!all(removed)) {
            stop("Failed to remove vignette HTML file(s):\n",
                 paste0("  ", html_files[!removed], collapse = "\n"),
                 call. = FALSE)
        }
    } else {
        message("No standalone HTML files found in vignettes/.")
    }

    pkgdown::build_site(pkg = pkg_root)

    llms_env <- new.env(parent = baseenv())
    source(file.path(pkg_root, "dev_scripts", "build_llms.R"),
           local = llms_env)
    llms_env$build_llms(pkg_root)

    invisible(pkg_root)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 1L) {
        stop("Usage: Rscript dev_scripts/build_site.R [package-root]",
             call. = FALSE)
    }
    build_mizer_site(if (length(args)) args[[1L]] else ".")
}
