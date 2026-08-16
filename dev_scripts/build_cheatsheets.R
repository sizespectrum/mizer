# Generate the cheatsheet vignettes and the upgrading article from the agent
# skills.
#
# Single source of truth: inst/skills/<topic>/SKILL.md. That file is shipped
# verbatim as a Claude Code / agent skill and is also the source for the
# corresponding article on the pkgdown website (the cheatsheets, plus
# vignettes/upgrading.Rmd). Never edit a generated vignette by hand -- edit the
# skill and re-run
#
#     source("dev_scripts/build_cheatsheets.R"); build_cheatsheets()
#
# What the generator does:
#   * strips the YAML frontmatter and the leading H1 (replaced by the vignette
#     title from `cheatsheet_topics` below)
#   * turns ```r into ```{r eval=FALSE}
#   * auto-links the first mention of each documented function to its pkgdown
#     reference page, using the alias -> .Rd mapping read from man/. In the
#     first column of a table the link is repeated on every row, since a table
#     is read as a lookup rather than in order.
#   * rewrites "the `foo` skill" cross-references into links to the matching
#     cheatsheet
#   * drops <!-- agent-only --> ... <!-- /agent-only --> blocks
#   * inserts a horizontal rule before each level-2 heading
#   * appends inst/skills/<topic>/quick-reference.md as a final section
#
# Code inside fenced blocks, inline code spans that are already links, and
# `$...$` math are never rewritten.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Topic table -----------------------------------------------------------------
# One entry per generated article, keyed by its vignette basename. `skill` names
# the directory under inst/skills/ whose SKILL.md is the article's source. The
# mapping is deliberately one-to-one: an article assembled from several skills
# cannot express a cross-reference to one of them, and collapses references to
# two different skills into the same link.
#
# Optional fields: `lead` is a paragraph inserted under the title; `setup` gives
# extra lines for the hidden setup chunk; `nolink` lists names never turned into
# reference links; `link_text` is how *other* articles refer to this one (it
# defaults to "<name> cheatsheet", which only reads correctly for a cheatsheet).
cheatsheet_topics <- list(
    `cheatsheet-model-setup` = list(
        skill  = "build-multispecies-model",
        title  = "Cheatsheet: Model Setup",
        nolink = "catchability",
        lead   = paste("This cheatsheet gives a quick overview of the functions",
                       "used to build a mizer model from a species parameter",
                       "data frame. For bringing the model to steady state and",
                       "calibrating it, see the",
                       "[calibration cheatsheet](cheatsheet-calibration.html)."),
        setup  = "params <- NS_params"
    ),
    `cheatsheet-calibration` = list(
        skill  = "calibrate-model",
        title  = "Cheatsheet: Steady State and Calibration",
        lead   = paste("This cheatsheet gives a quick overview of the functions",
                       "used to bring a mizer model to steady state and calibrate",
                       "it to observed biomasses, yields and growth.",
                       "For full documentation of each function, follow the links."),
        setup  = "params <- NS_params"
    ),
    `cheatsheet-changing-parameters` = list(
        skill  = "change-parameters",
        title  = "Cheatsheet: Changing Model Parameters",
        lead   = "",
        setup  = "params <- NS_params",
        nolink = c("selectivity", "catchability", "maturity")
    ),
    `cheatsheet-fishing` = list(
        skill  = "set-up-fishing",
        title  = "Cheatsheet: Fishing",
        lead   = paste("This cheatsheet gives a quick overview of how fishing",
                       "is set up in mizer: gears, selectivity, catchability,",
                       "and effort. For full documentation of each function,",
                       "follow the links."),
        setup  = c("params <- NS_params", "sim <- NS_sim"),
        # Names that double as gear_params/species_params column names: linking
        # the bare mention would point at a function the reader did not mean.
        nolink = c("catchability", "selectivity")
    ),
    `cheatsheet-running-simulations` = list(
        skill  = "run-simulation",
        title  = "Cheatsheet: Running Simulations",
        lead   = paste("This cheatsheet covers projecting a model forward with",
                       "project(): its arguments, the four ways of specifying",
                       "fishing effort, continuing and comparing runs, and the",
                       "choice of numerical scheme."),
        setup  = c("params <- NS_params", "sim <- NS_sim")
    ),
    `cheatsheet-analysis-and-plotting` = list(
        skill  = "analyse-and-plot",
        title  = "Cheatsheet: Analysis and Plotting",
        # Base generics: a reader clicking these expects base R, not mizer.
        nolink = c("print", "as.data.frame"),
        lead   = paste("This cheatsheet gives a quick overview of the functions",
                       "available in mizer for analysing the results of",
                       "simulations and creating plots. For full documentation",
                       "of each function, follow the links."),
        setup  = c("params <- NS_params", "sim <- NS_sim")
    ),
    `cheatsheet-stability` = list(
        skill  = "analyse-stability",
        title  = "Cheatsheet: Dynamic Stability",
        lead   = paste("This cheatsheet gives a quick overview of the",
                       "experimental tools for analysing the dynamic stability",
                       "of a mizer steady state and characterising the limit",
                       "cycles that can replace it. For full documentation of",
                       "each function, follow the links."),
        setup  = "params <- NS_params"
    ),
    `cheatsheet-extending-mizer` = list(
        skill  = "extend-mizer",
        title  = "Cheatsheet: Extending mizer",
        lead   = paste("This cheatsheet gives a quick overview of the mechanisms",
                       "for customising mizer's dynamics, from external",
                       "encounter and mortality through to whole new ecosystem",
                       "components. For the full treatment see the",
                       "[Extending mizer](extending-mizer.html) article."),
        setup  = "params <- NS_params"
    ),
    # Not a cheatsheet: the upgrading article is the release-by-release list of
    # changes that break existing code, and the skill an agent loads when a
    # user's script stops working after an upgrade. The skill's symptom index,
    # which is of no use to a human reading by release, is agent-only.
    upgrading = list(
        skill     = "upgrade-mizer-code",
        title     = "Upgrading mizer",
        link_text = "Upgrading mizer article",
        # Deliberately no `lead` and no `setup`: the skill body opens with its
        # own introduction, and nothing in the article is evaluated.
        nolink    = c("print", "summary", "plot", "as.data.frame")
    )
)

# Map a skill name to the article it becomes, and to the words other articles
# use to refer to it, for cross-references.
skill_links <- local({
    out <- list()
    for (vig in names(cheatsheet_topics)) {
        spec <- cheatsheet_topics[[vig]]
        out[[spec$skill]] <- list(
            vignette = vig,
            text = spec$link_text %||%
                paste(sub("^cheatsheet-", "", vig), "cheatsheet")
        )
    }
    out
})


# Alias -> reference page ------------------------------------------------------

#' Build a map from documented alias to pkgdown reference page basename
#'
#' pkgdown writes one `reference/<Rd basename>.html` per `.Rd` file, reachable
#' under every alias that file declares. So `catchability` maps to `setFishing`.
rd_alias_map <- function(man_dir = "man") {
    rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
    map <- new.env(parent = emptyenv())
    for (f in rd_files) {
        topic <- sub("\\.Rd$", "", basename(f))
        lines <- readLines(f, warn = FALSE)
        aliases <- sub("^\\\\alias\\{(.*)\\}\\s*$", "\\1", grep("^\\\\alias\\{", lines, value = TRUE))
        for (a in aliases) {
            # An alias may be escaped in the Rd, e.g. \alias{\%*\%}
            a <- gsub("\\\\", "", a)
            if (!nzchar(a)) next
            # First .Rd wins, but prefer the file named after the alias itself.
            if (is.null(map[[a]]) || identical(topic, a)) map[[a]] <- topic
        }
    }
    map
}


# Inline-code linking ----------------------------------------------------------

#' Wrap the first mention of each documented function in a pkgdown link
#'
#' Only inline code spans are considered. A span is linked when it is a call
#' such as `getFMort(params)`, or a bare identifier of at least `min_bare`
#' characters such as `gear_params`. Short bare names (`w`, `N`) read as prose
#' rather than as cross-references and are skipped, as are the bare names listed
#' in `nolink`, which double as data-frame column names.
#'
#' `force_to` suspends the first-mention rule up to that character position: a
#' name is linked there even if it was linked earlier in the article. It is used
#' for the first column of a table, which is read as a lookup, out of order and
#' without the prose that introduced the name.
link_first_mentions <- function(line, map, seen, nolink = character(0),
                                min_bare = 3L, force_to = 0L) {
    m <- gregexpr("`[^`]+`", line, perl = TRUE)[[1]]
    if (m[1] == -1) return(list(line = line, seen = seen))
    starts <- as.integer(m)
    lens <- attr(m, "match.length")

    # Pass 1, left to right: decide which spans to link, so that the *first*
    # mention on a line wins rather than the last.
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
            # Bare identifier: must stand alone and be long enough to read as a
            # cross-reference rather than as prose.
            if (!identical(name, inner)) next          # e.g. `gp$l50`
            if (nchar(name) < min_bare) next           # e.g. `w`, `N`
        }

        if (is.null(map[[name]])) next
        if (isTRUE(seen[[name]]) && en > force_to) next
        seen[[name]] <- TRUE
        todo <- c(todo, i)
    }
    if (!length(todo)) return(list(line = line, seen = seen))

    # Pass 2, right to left, so the earlier match positions stay valid.
    out <- line
    for (i in rev(todo)) {
        st <- starts[i]
        en <- st + lens[i] - 1L
        span <- substr(out, st, en)
        name <- sub("^`([A-Za-z._][A-Za-z0-9._]*).*$", "\\1", span)
        repl <- sprintf("[%s](../reference/%s.html)", span, map[[name]])
        out <- paste0(substr(out, 1L, st - 1L), repl,
                      if (en < nchar(out)) substr(out, en + 1L, nchar(out)) else "")
    }
    list(line = out, seen = seen)
}


#' Character position up to which linking is forced (all table rows)
#'
#' Returns 0 for anything that is not a markdown table row, and nchar(line) for
#' a table row, so all cells in a table are always linked. The delimiter row is
#' excluded.
table_key_cell <- function(line) {
    if (!grepl("^\\s*\\|", line)) return(0L)
    if (grepl("^\\s*\\|[-: |]+\\|\\s*$", line)) return(0L)
    nchar(line)
}


# Main conversion --------------------------------------------------------------

#' Read one SKILL.md and return its body: frontmatter, H1 and agent-only
#' blocks removed.
skill_body <- function(skill, pkg_root = ".") {
    lines <- readLines(file.path(pkg_root, "inst", "skills", skill, "SKILL.md"),
                       warn = FALSE)

    # Strip YAML frontmatter.
    if (length(lines) && grepl("^---\\s*$", lines[1])) {
        close <- which(grepl("^---\\s*$", lines))[2]
        lines <- lines[-seq_len(close)]
    }
    # Strip the leading H1 (the vignette title replaces it) and any blank lines
    # left in front of the body.
    first <- which(nzchar(trimws(lines)))[1]
    if (!is.na(first) && grepl("^# ", lines[first])) lines <- lines[-seq_len(first)]
    first <- which(nzchar(trimws(lines)))[1]
    if (!is.na(first) && first > 1L) lines <- lines[-seq_len(first - 1L)]
    lines
}

skill_to_cheatsheet <- function(vignette, spec, map, pkg_root = ".") {
    lines <- skill_body(spec$skill, pkg_root)

    # Drop agent-only blocks.
    keep <- rep(TRUE, length(lines))
    inside <- FALSE
    for (i in seq_along(lines)) {
        if (grepl("<!--\\s*agent-only\\s*-->", lines[i])) inside <- TRUE
        if (inside) keep[i] <- FALSE
        if (grepl("<!--\\s*/agent-only\\s*-->", lines[i])) inside <- FALSE
    }
    lines <- lines[keep]

    # Walk the body: convert fences, link prose, add rules before H2s.
    out <- character(0)
    in_fence <- FALSE
    seen <- list()
    for (ln in lines) {
        if (grepl("^\\s*```", ln)) {
            if (!in_fence) {
                ln <- sub("^(\\s*```)r\\s*$", "\\1{r eval=FALSE}", ln)
                in_fence <- TRUE
            } else {
                in_fence <- FALSE
            }
            out <- c(out, ln)
            next
        }
        if (in_fence) {
            out <- c(out, ln)
            next
        }
        if (grepl("^## ", ln)) out <- c(out, "---", "")
        # Headings are never linked: a link inside a heading nests inside the
        # table-of-contents entry, which is itself a link. The name keeps its
        # first-mention link, which lands on the prose below instead.
        if (grepl("^#{1,6} ", ln)) {
            out <- c(out, ln)
            next
        }
        # Cross-references to sibling skills become links to their articles.
        for (s in names(skill_links)) {
            target <- skill_links[[s]]
            ln <- gsub(sprintf("`%s` skill", s),
                       sprintf("[%s](%s.html)", target$text, target$vignette),
                       ln, fixed = TRUE)
        }
        res <- link_first_mentions(ln, map, seen,
                                   nolink = spec$nolink %||% character(0),
                                   force_to = table_key_cell(ln))
        seen <- res$seen
        out <- c(out, res$line)
    }

    # Append the quick reference section, if the skill ships one.
    qr_path <- file.path(pkg_root, "inst", "skills", spec$skill,
                         "quick-reference.md")
    if (file.exists(qr_path)) {
        qr <- readLines(qr_path, warn = FALSE)
        qr <- sub("^(\\s*```)r\\s*$", "\\1{r eval=FALSE}", qr)
        out <- c(out, "", "---", "", "## Quick reference", "", qr)
    }

    sources <- paste0("inst/skills/", spec$skill, "/SKILL.md")
    header <- c(
        "---",
        sprintf("title: \"%s\"", spec$title),
        "output:",
        "  html_document:",
        "    toc: yes",
        "    toc_float: yes",
        "    fig_width: 5",
        "    fig_height: 4",
        "---",
        "",
        paste0("<!-- Generated from ", sources,
               " by dev_scripts/build_cheatsheets.R."),
        "     Do not edit this file by hand -- edit the skill and re-run the generator. -->",
        "",
        "```{r include=FALSE}",
        "library(mizer)",
        spec$setup,
        "```",
        ""
    )
    if (!is.null(spec$lead)) {
        header <- c(header, strwrap(spec$lead, width = 78), "")
    }

    c(header, out)
}

#' Warn about links that will 404 and about skill cross-references the rewriter
#' did not convert
#'
#' Three ways a generated article can point at nothing: a `../reference/` link
#' with no matching man page; a sibling-article link with no matching vignette
#' (these come from hand-written `lead` text and from the skills' own article
#' links, so nothing else validates them); and a leftover `` `<skill>` `` code
#' span, meaning the source phrased the cross-reference in some way other than
#' "`<skill>` skill" and a reader is being sent to a skill they cannot see.
check_output <- function(txt, vignette, man_dir, vig_dir) {
    pat <- "\\.\\./reference/[A-Za-z0-9_.]+\\.html"
    hits <- unlist(regmatches(txt, gregexpr(pat, txt)))
    targets <- unique(sub("\\.html$", "",
                          sub("^\\.\\./reference/", "", hits)))
    dead <- targets[!file.exists(file.path(man_dir, paste0(targets, ".Rd")))]
    if (length(dead)) {
        warning(vignette, ": link target(s) with no man page: ",
                paste(dead, collapse = ", "), call. = FALSE)
    }

    # Links to sibling articles: "](name.html)" with no directory part.
    apat <- "\\]\\([A-Za-z0-9_-]+\\.html(#[A-Za-z0-9_-]+)?\\)"
    ahits <- unlist(regmatches(txt, gregexpr(apat, txt)))
    arts <- unique(sub("\\.html.*$", "", sub("^\\]\\(", "", ahits)))
    exists_art <- vapply(arts, function(a) {
        any(file.exists(file.path(vig_dir, paste0(a, c(".Rmd", ".qmd", ".rmd")))))
    }, logical(1))
    dead_art <- arts[!exists_art]
    if (length(dead_art)) {
        warning(vignette, ": article link(s) with no vignette: ",
                paste(dead_art, collapse = ", "), call. = FALSE)
        dead <- c(dead, dead_art)
    }

    body <- paste(txt, collapse = "\n")
    is_stray <- vapply(names(skill_links), function(s) {
        grepl(paste0("`", s, "`"), body, fixed = TRUE)
    }, logical(1))
    stray <- names(skill_links)[is_stray]
    if (length(stray)) {
        warning(vignette, ": unconverted skill reference(s): ",
                paste(stray, collapse = ", "),
                ". Phrase each as \"`<skill-name>` skill\" in the source.",
                call. = FALSE)
    }
    length(dead) + length(stray)
}

build_cheatsheets <- function(pkg_root = ".",
                              vignettes = names(cheatsheet_topics)) {
    man_dir <- file.path(pkg_root, "man")
    vig_dir <- file.path(pkg_root, "vignettes")
    map <- rd_alias_map(man_dir)
    bad <- 0L
    for (vig in vignettes) {
        spec <- cheatsheet_topics[[vig]]
        txt <- skill_to_cheatsheet(vig, spec, map, pkg_root)
        bad <- bad + check_output(txt, vig, man_dir, vig_dir)
        dest <- file.path(pkg_root, "vignettes", paste0(vig, ".Rmd"))
        writeLines(txt, dest)
        message("Wrote ", dest, " from ", spec$skill)
    }
    if (bad == 0L) {
        message("All links resolve, all skill references converted.")
    }
    invisible(TRUE)
}
