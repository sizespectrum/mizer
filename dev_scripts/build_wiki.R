# Generate GitHub wiki documentation from the developer agent skills in .claude/skills/
#
# Single source of truth: .claude/skills/<topic>.md.
# These files serve as instructions for AI agents and are converted into human
# developer guides on the mizer GitHub wiki (in ../mizer.wiki/).
#
# Usage:
#     source("dev_scripts/build_wiki.R"); build_wiki()
#
# Or from shell:
#     Rscript dev_scripts/build_wiki.R
#     Rscript dev_scripts/build_wiki.R --wiki ../mizer.wiki

`%||%` <- function(x, y) if (is.null(x)) y else x

# Topic table -----------------------------------------------------------------
# One entry per generated wiki page, keyed by its wiki page name. `skill` names
# the file in .claude/skills/<skill>.md.
wiki_topics <- list(
    `Array-Wrapper-Classes` = list(
        skill = "array-wrapper-classes",
        title = "Array Wrapper Classes",
        nolink = c("ArraySpeciesBySize", "ArrayTimeBySpecies",
                   "ArrayTimeBySpeciesBySize", "ArrayResourceBySize",
                   "ArrayTimeByResourceBySize")
    ),
    `Size-Grid-Integrals` = list(
        skill = "size-grid-integrals",
        title = "Integrating over the Size Grid",
        nolink = c("w", "N", "dw")
    ),
    `Species-Parameter-Defaults` = list(
        skill = "species-param-defaults",
        title = "Species Parameter Defaults",
        nolink = c("catchability", "selectivity", "maturity")
    ),
    `Condition-and-Signal-Handling` = list(
        skill = "info-signals",
        title = "Condition and Signal Handling"
    ),
    `Documenting-S3-Generics` = list(
        skill = "document-s3-generics",
        title = "Documenting S3 Generics",
        nolink = c("print", "summary", "plot", "as.data.frame", "str")
    ),
    `Upgrading-Mizer-Data-Objects` = list(
        skill = "upgrade-mizer-data",
        title = "Upgrading Mizer Data Objects"
    ),
    `Test-Organisation` = list(
        skill = "test-organisation",
        title = "Test Organisation"
    ),
    `Test-Fixtures` = list(
        skill = "test-fixtures",
        title = "Test Fixtures and Keeping the Suite Fast"
    ),
    `Building-Documentation-and-Website` = list(
        skill = "build-documentation",
        title = "Building Documentation and Website"
    )
)

# User skill links (on pkgdown site)
user_skill_links <- list(
    `build-multispecies-model` = list(
        text = "Model Setup cheatsheet",
        url = "https://sizespectrum.org/mizer/articles/cheatsheet-model-setup.html"
    ),
    `calibrate-model` = list(
        text = "Calibration cheatsheet",
        url = "https://sizespectrum.org/mizer/articles/cheatsheet-calibration.html"
    ),
    `change-parameters` = list(
        text = "Changing Parameters cheatsheet",
        url = "https://sizespectrum.org/mizer/articles/cheatsheet-changing-parameters.html"
    ),
    `set-up-fishing` = list(
        text = "Fishing cheatsheet",
        url = "https://sizespectrum.org/mizer/articles/cheatsheet-fishing.html"
    ),
    `run-simulation` = list(
        text = "Running Simulations cheatsheet",
        url = "https://sizespectrum.org/mizer/articles/cheatsheet-running-simulations.html"
    ),
    `analyse-and-plot` = list(
        text = "Analysis and Plotting cheatsheet",
        url = "https://sizespectrum.org/mizer/articles/cheatsheet-analysis-and-plotting.html"
    ),
    `analyse-stability` = list(
        text = "Dynamic Stability cheatsheet",
        url = "https://sizespectrum.org/mizer/articles/cheatsheet-stability.html"
    ),
    `extend-mizer` = list(
        text = "Extending mizer cheatsheet",
        url = "https://sizespectrum.org/mizer/articles/cheatsheet-extending-mizer.html"
    ),
    `upgrade-mizer-code` = list(
        text = "Upgrading mizer guide",
        url = "https://sizespectrum.org/mizer/articles/upgrading.html"
    )
)

# Alias -> reference page ------------------------------------------------------

rd_alias_map <- function(man_dir = "man") {
    rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
    map <- new.env(parent = emptyenv())
    for (f in rd_files) {
        topic <- sub("\\.Rd$", "", basename(f))
        lines <- readLines(f, warn = FALSE)
        aliases <- sub("^\\\\alias\\{(.*)\\}\\s*$", "\\1", grep("^\\\\alias\\{", lines, value = TRUE))
        for (a in aliases) {
            a <- gsub("\\\\", "", a)
            if (!nzchar(a)) next
            if (is.null(map[[a]]) || identical(topic, a)) map[[a]] <- topic
        }
    }
    map
}

# Inline-code linking ----------------------------------------------------------

link_first_mentions <- function(line, map, seen, nolink = character(0),
                                min_bare = 3L, force_to = 0L,
                                url_prefix = "https://sizespectrum.org/mizer/reference/") {
    m <- gregexpr("`[^`]+`", line, perl = TRUE)[[1]]
    if (m[1] == -1) return(list(line = line, seen = seen))
    starts <- as.integer(m)
    lens <- attr(m, "match.length")

    todo <- integer(0)
    for (i in seq_along(starts)) {
        st <- starts[i]
        en <- st + lens[i] - 1L
        inner <- substr(line, st + 1L, en - 1L)

        # Already the text of a markdown link: "[`foo`](...)"
        before <- if (st > 1L) substr(line, st - 1L, st - 1L) else ""
        after <- if (en < nchar(line)) substr(line, en + 1L, en + 1L) else ""
        if (identical(before, "[") && identical(after, "]")) next

        name <- sub("^([A-Za-z._][A-Za-z0-9._]*).*$", "\\1", inner)
        if (!nzchar(name) || !grepl("^[A-Za-z._]", inner)) next
        if (name %in% nolink) next

        if (!grepl("^[A-Za-z._][A-Za-z0-9._]*\\(", inner)) {
            if (!identical(name, inner)) next
            if (nchar(name) < min_bare) next
        }

        if (is.null(map[[name]])) next
        if (isTRUE(seen[[name]]) && en > force_to) next
        seen[[name]] <- TRUE
        todo <- c(todo, i)
    }
    if (!length(todo)) return(list(line = line, seen = seen))

    out <- line
    for (i in rev(todo)) {
        st <- starts[i]
        en <- st + lens[i] - 1L
        span <- substr(out, st, en)
        name <- sub("^`([A-Za-z._][A-Za-z0-9._]*).*$", "\\1", span)
        repl <- sprintf("[%s](%s%s.html)", span, url_prefix, map[[name]])
        out <- paste0(substr(out, 1L, st - 1L), repl,
                      if (en < nchar(out)) substr(out, en + 1L, nchar(out)) else "")
    }
    list(line = out, seen = seen)
}

table_key_cell <- function(line) {
    if (!grepl("^\\s*\\|", line)) return(0L)
    if (grepl("^\\s*\\|[-: |]+\\|\\s*$", line)) return(0L)
    close <- regexpr("\\|[^|]*\\|", line, perl = TRUE)
    if (close == -1L) return(0L)
    close + attr(close, "match.length") - 1L
}

# Conversion logic -------------------------------------------------------------

skill_to_wiki <- function(page_name, spec, map, pkg_root = ".") {
    skill_file <- file.path(pkg_root, ".claude", "skills", paste0(spec$skill, ".md"))
    if (!file.exists(skill_file)) {
        stop("Skill file not found: ", skill_file)
    }
    lines <- readLines(skill_file, warn = FALSE)

    # Drop agent-only blocks
    keep <- rep(TRUE, length(lines))
    inside <- FALSE
    for (i in seq_along(lines)) {
        if (grepl("<!--\\s*agent-only\\s*-->", lines[i])) inside <- TRUE
        if (inside) keep[i] <- FALSE
        if (grepl("<!--\\s*/agent-only\\s*-->", lines[i])) inside <- FALSE
    }
    lines <- lines[keep]

    out <- character(0)
    in_fence <- FALSE
    seen <- list()

    for (ln in lines) {
        if (grepl("^\\s*```", ln)) {
            in_fence <- !in_fence
            out <- c(out, ln)
            next
        }
        if (in_fence) {
            out <- c(out, ln)
            next
        }

        # Headings: don't link inside headings
        if (grepl("^#{1,6} ", ln)) {
            out <- c(out, ln)
            next
        }

        # Rewrite developer skill cross-references (.claude/skills/<topic>.md or `<topic>` skill)
        for (w_name in names(wiki_topics)) {
            w_spec <- wiki_topics[[w_name]]
            # Pattern 1: `.claude/skills/foo.md`
            ln <- gsub(paste0("`\\.claude/skills/", w_spec$skill, "\\.md`"),
                       sprintf("[%s](%s)", w_spec$title, w_name),
                       ln)
            ln <- gsub(paste0("\\.claude/skills/", w_spec$skill, "\\.md"),
                       sprintf("[%s](%s)", w_spec$title, w_name),
                       ln)
            # Pattern 2: `foo` skill
            ln <- gsub(sprintf("`%s` skill", w_spec$skill),
                       sprintf("[%s](%s)", w_spec$title, w_name),
                       ln, fixed = TRUE)
        }

        # Rewrite user skill cross-references
        for (u_name in names(user_skill_links)) {
            u_spec <- user_skill_links[[u_name]]
            ln <- gsub(sprintf("`%s` skill", u_name),
                       sprintf("[%s](%s)", u_spec$text, u_spec$url),
                       ln, fixed = TRUE)
        }

        # Rewrite relative article links like `[Extending mizer](extending-mizer.html)` or `vignettes/foo.qmd`
        ln <- gsub("\\]\\(([a-zA-Z0-9_-]+)\\.html(#[a-zA-Z0-9_-]+)?\\)",
                   "](https://sizespectrum.org/mizer/articles/\\1.html\\2)",
                   ln)
        ln <- gsub("`vignettes/([a-zA-Z0-9_-]+)\\.(qmd|Rmd)`",
                   "[`\\1`](https://sizespectrum.org/mizer/articles/\\1.html)",
                   ln)

        res <- link_first_mentions(ln, map, seen,
                                   nolink = spec$nolink %||% character(0),
                                   force_to = table_key_cell(ln),
                                   url_prefix = "https://sizespectrum.org/mizer/reference/")
        seen <- res$seen
        out <- c(out, res$line)
    }

    banner <- c(
        paste0("<!-- Generated from .claude/skills/", spec$skill,
               ".md by dev_scripts/build_wiki.R. -->"),
        "<!-- Do not edit this file by hand -- edit the skill in mizer and re-run the generator. -->",
        ""
    )

    c(banner, out)
}

generate_sidebar <- function() {
    c(
        "<!-- Generated by dev_scripts/build_wiki.R -->",
        "### mizer Developer Wiki",
        "",
        "- [[Home]]",
        "- [Developer Guide (Overview)](Developer-Guide)",
        "",
        "#### Contributing & Setup",
        "- [[Working with git and GitHub|Working-with-git-and-GitHub]]",
        "- [[Developer FAQ|Developer-FAQ]]",
        "",
        "#### Architecture & Technical Guides",
        "- [[Array Wrapper Classes|Array-Wrapper-Classes]]",
        "- [[Integrating over Size Grid|Size-Grid-Integrals]]",
        "- [[Species Parameter Defaults|Species-Parameter-Defaults]]",
        "- [[Condition & Signal Handling|Condition-and-Signal-Handling]]",
        "- [[Documenting S3 Generics|Documenting-S3-Generics]]",
        "- [[Upgrading Mizer Data Objects|Upgrading-Mizer-Data-Objects]]",
        "- [[Test Organisation|Test-Organisation]]",
        "- [[Test Fixtures & Speed|Test-Fixtures]]",
        "- [[Building Docs & Website|Building-Documentation-and-Website]]",
        "",
        "---",
        "**Online Links**",
        "- [mizer Website](https://sizespectrum.org/mizer/)",
        "- [GitHub Repository](https://github.com/sizespectrum/mizer)"
    )
}

generate_home <- function() {
    c(
        "# Welcome to the mizer Developer Wiki",
        "",
        "This wiki contains technical architecture documentation and guides for developers, maintainers, and contributors to the [mizer](https://github.com/sizespectrum/mizer) R package.",
        "",
        "For user guides, cheatsheets, and function references, visit the [mizer website](https://sizespectrum.org/mizer/).",
        "",
        "---",
        "",
        "## Getting Started as a Contributor",
        "",
        "- **[Working with git and GitHub](Working-with-git-and-GitHub)**: Step-by-step tutorial on branches, commits, pull requests, and syncing with upstream.",
        "- **[Developer FAQ](Developer-FAQ)**: Common contributor questions (branch management, clean master branch, etc.).",
        "",
        "---",
        "",
        "## Architecture & Developer Guides",
        "",
        "- **[Array Wrapper Classes](Array-Wrapper-Classes)**: The `ArraySpeciesBySize` and sibling S3 wrapper classes, `Ops` arithmetic stripping, and slot assignment rules.",
        "- **[Integrating over the Size Grid](Size-Grid-Integrals)**: First-order vs second-order quadrature schemes (`second_order_w`), `sizeIntegral()`, and avoiding double-counting traps.",
        "- **[Species Parameter Defaults](Species-Parameter-Defaults)**: Where parameter defaults belong (setter-owned vs central), and `given_species_params` mechanics.",
        "- **[Condition and Signal Handling](Condition-and-Signal-Handling)**: The `signal_info()` system, `info_level`, severity levels, and frozen rate warnings.",
        "- **[Documenting S3 Generics](Documenting-S3-Generics)**: Documenting multi-method S3 generics on a shared man page without `R CMD check` codoc warnings.",
        "- **[Upgrading Mizer Data Objects](Upgrading-Mizer-Data-Objects)**: Steps for adding slots to `MizerParams` / `MizerSim`, updating `upgradeParams()`, and bumping version gating.",
        "- **[Test Organisation](Test-Organisation)**: Test file naming conventions (`test-<file>.R`) and snapshot rules.",
        "- **[Test Fixtures and Keeping the Suite Fast](Test-Fixtures)**: Shared `_small` fixtures in `helper.R`, `delayedAssign()`, and parallel test execution.",
        "- **[Building Documentation and Website](Building-Documentation-and-Website)**: Generating cheatsheets, maintaining `llms.txt`, and configuring `_pkgdown.yml`."
    )
}

generate_developer_guide_overview <- function() {
    c(
        "# mizer Developer Guide",
        "",
        "This guide provides an overview of mizer's internal architecture, design principles, and developer workflows.",
        "",
        "## Architecture Overview",
        "",
        "- **`MizerParams` (S4)**: Central object representing the model specification. Modified via setter functions (`setFishing()`, `setMetabolicRate()`, etc.) that return new copies.",
        "- **`MizerSim` (S4)**: Holds simulation time-series results (`n`, `n_pp`, `n_other`, `effort`).",
        "- **S3 Methods on S4 Classes**: Mizer deliberately registers methods for `MizerParams` and `MizerSim` as S3 methods (e.g. `plot.MizerParams`) rather than `setMethod()`, to avoid promoting base generics to S4 generics session-wide.",
        "- **Two Quadrature Schemes**: `second_order_w` selects first-order boundary values vs second-order finite-volume bin averages.",
        "- **Gated Validation**: `validParams()` uses `validation_key()` fingerprinting to cache validation checks efficiently.",
        "",
        "---",
        "",
        "## Detailed Developer Guides",
        "",
        "| Guide | Description |",
        "|---|---|",
        "| **[Array Wrapper Classes](Array-Wrapper-Classes)** | S3 wrapper classes around rate arrays, `slot[] <- value` rules, and method outputs |",
        "| **[Integrating over Size Grid](Size-Grid-Integrals)** | Bin averaging, `sizeIntegral()`, encounter kernels, and quadrature rules |",
        "| **[Species Parameter Defaults](Species-Parameter-Defaults)** | Ownership of parameter defaults (setter-owned vs central) and `given_species_params` |",
        "| **[Condition & Signal Handling](Condition-and-Signal-Handling)** | User notifications, `info_level`, `signal_info()`, and frozen parameter warnings |",
        "| **[Documenting S3 Generics](Documenting-S3-Generics)** | Shared man pages, `@usage NULL`, `@param ...`, and avoiding `codoc` check errors |",
        "| **[Upgrading Mizer Data](Upgrading-Mizer-Data-Objects)** | Class changes, new slots, updating `upgradeParams()`, and updating stored fixtures |",
        "| **[Test Organisation](Test-Organisation)** | Naming conventions for test files, snapshots, and splitting files |",
        "| **[Test Fixtures & Speed](Test-Fixtures)** | Shared `_small` fixtures, parallel workers, and experimental test gating |",
        "| **[Building Documentation & Website](Building-Documentation-and-Website)** | Cheatsheet generation, `llms.txt`, `_pkgdown.yml`, and navbar syncing |"
    )
}

# Main build function ----------------------------------------------------------

build_wiki <- function(pkg_root = ".", wiki_root = "../mizer.wiki") {
    man_dir <- file.path(pkg_root, "man")
    if (!dir.exists(wiki_root)) {
        stop("Wiki directory not found at: ", wiki_root)
    }

    map <- rd_alias_map(man_dir)
    message("Building wiki pages from .claude/skills/ into ", wiki_root, "...")

    for (page_name in names(wiki_topics)) {
        spec <- wiki_topics[[page_name]]
        txt <- skill_to_wiki(page_name, spec, map, pkg_root)
        dest <- file.path(wiki_root, paste0(page_name, ".md"))
        writeLines(txt, dest)
        message("  Wrote ", basename(dest), " (from .claude/skills/", spec$skill, ".md)")
    }

    # Write _Sidebar.md
    sidebar_dest <- file.path(wiki_root, "_Sidebar.md")
    writeLines(generate_sidebar(), sidebar_dest)
    message("  Wrote _Sidebar.md")

    # Write Home.md
    home_dest <- file.path(wiki_root, "Home.md")
    writeLines(generate_home(), home_dest)
    message("  Wrote Home.md")

    # Write Developer-Guide.md
    devguide_dest <- file.path(wiki_root, "Developer-Guide.md")
    writeLines(generate_developer_guide_overview(), devguide_dest)
    message("  Wrote Developer-Guide.md")

    message("Wiki build complete.")
    invisible(TRUE)
}

# Allow direct CLI invocation
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    wiki_dir <- "../mizer.wiki"
    if (length(args) >= 2 && args[1] == "--wiki") {
        wiki_dir <- args[2]
    }
    build_wiki(wiki_root = wiki_dir)
}
