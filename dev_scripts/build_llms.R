# Finish the llms.txt file that pkgdown generates, and ship it with the package.
#
# Run from the mizer repo root after every pkgdown::build_site():
#
#     source("dev_scripts/build_llms.R"); build_llms()
#
# pkgdown writes docs/llms.txt by prepending a rendering of README.md to the
# categorised package index. That preamble is full of badge and logo markdown,
# which is noise to an agent and is not what the llms.txt convention asks for --
# a short prose orientation followed by the links. This script replaces it with
# the hand-written header in pkgdown/llms-header.md, which is the file to edit
# when the orientation needs to change.
#
# The result is written to two places:
#
#   * docs/llms.txt   -- served at https://sizespectrum.org/mizer/llms.txt
#   * inst/llms.txt   -- installed with the package, so that
#                        mizerAgents::setup_mizer_agent() can point an agent at
#                        the index for the mizer it actually has installed.
#                        (docs/ is .Rbuildignore'd, so the copy is not optional.)
#
# The script is idempotent: docs/llms.txt is rewritten in place, so it must be
# safe to run twice. Everything above the "# Package index" heading is replaced,
# and that heading is what both the freshly generated and the already-processed
# file are cut at.

INDEX_HEADING <- "# Package index"

build_llms <- function(pkg_root = ".") {
    src    <- file.path(pkg_root, "docs", "llms.txt")
    header <- file.path(pkg_root, "pkgdown", "llms-header.md")

    if (!file.exists(src)) {
        stop("`", src, "` not found. Run pkgdown::build_site() first.",
             call. = FALSE)
    }
    if (!file.exists(header)) {
        stop("`", header, "` not found.", call. = FALSE)
    }

    lines <- readLines(src, warn = FALSE)
    start <- which(lines == INDEX_HEADING)
    if (length(start) != 1L) {
        stop("Expected exactly one \"", INDEX_HEADING, "\" line in ", src,
             ", found ", length(start),
             ". The pkgdown output has changed shape; check it by hand.",
             call. = FALSE)
    }

    out <- c(readLines(header, warn = FALSE), "", lines[start:length(lines)])

    for (dest in c(src, file.path(pkg_root, "inst", "llms.txt"))) {
        writeLines(out, dest)
        message("Wrote ", dest, " (", length(out), " lines)")
    }

    invisible(out)
}
