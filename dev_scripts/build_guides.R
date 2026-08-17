# Generate the guide vignettes and the upgrading article from the agent
# skills.
#
# Single source of truth: inst/skills/<topic>/SKILL.md. That file is shipped
# verbatim as a Claude Code / agent skill and is also the source for the
# corresponding article on the pkgdown website (the guides, plus
# vignettes/upgrading.Rmd). Never edit a generated vignette by hand -- edit the
# skill and re-run
#
#     source("dev_scripts/build_guides.R"); build_guides()
#
# What the generator does:
#   * strips the YAML frontmatter and the leading H1 (replaced by the vignette
#     title from `guide_topics` below)
#   * turns ```r into ```{r eval=FALSE}
#   * auto-links the first mention of each documented function to its pkgdown
#     reference page, using the alias -> .Rd mapping read from man/. In the
#     first column of a table the link is repeated on every row, since a table
#     is read as a lookup rather than in order.
#   * rewrites "the `foo` skill" cross-references into links to the matching
#     guide
#   * drops <!-- agent-only --> ... <!-- /agent-only --> blocks
#   * inserts a horizontal rule before each level-2 heading
#   * appends inst/skills/<topic>/quick-reference.md as a final section
#
# Code inside fenced blocks, inline code spans that are already links, and
# `$...$` math are never rewritten.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Topic table -----------------------------------------------------------------
# One entry per generated article, keyed by the directory under inst/skills/
# whose SKILL.md is its source. That skill name is the one name for the topic:
# the article is `guide-<skill>.Rmd` and its title is "Guide: " followed by the
# skill's own H1, so neither can drift from the skill. The mapping is
# deliberately one-to-one: an article assembled from several skills cannot
# express a cross-reference to one of them, and collapses references to two
# different skills into the same link.
#
# Optional fields: `lead` is a paragraph inserted under the title; `setup` gives
# extra lines for the hidden setup chunk; `nolink` lists names never turned into
# reference links. `vignette`, `title` and `link_text` override the derived
# values and should be needed only for an article that is not a guide.
guide_topics <- list(
    `understand-size-spectrum-dynamics` = list(
        lead   = paste("This guide gives a concise overview of the core",
                       "principles of size-spectrum modelling in mizer: physiology,",
                       "emergent growth and mortality, density dependence, and",
                       "trophic feedbacks. For the full mathematical formulation,",
                       "see the [model description](model_description.html) article."),
        setup  = "params <- NS_params"
    ),
    `build-model` = list(
        nolink = "catchability",
        lead   = paste("This guide gives an overview of the functions",
                       "used to build a mizer model from a species parameter",
                       "data frame. For bringing the model to steady state and",
                       "calibrating it, see the",
                       "[guide to reaching steady state and calibrating](guide-calibrate-model.html)."),
        setup  = "params <- NS_params"
    ),
    `calibrate-model` = list(
        lead   = paste("This guide gives an overview of the functions",
                       "used to bring a mizer model to steady state and calibrate",
                       "it to observed biomasses, yields and growth.",
                       "For full documentation of each function, follow the links."),
        setup  = "params <- NS_params"
    ),
    `change-parameters` = list(
        lead   = "",
        setup  = "params <- NS_params",
        nolink = c("selectivity", "catchability", "maturity")
    ),
    `set-up-fishing` = list(
        lead   = paste("This guide gives an overview of how fishing",
                       "is set up in mizer: gears, selectivity, catchability,",
                       "and effort. For full documentation of each function,",
                       "follow the links."),
        setup  = c("params <- NS_params", "sim <- NS_sim"),
        # Names that double as gear_params/species_params column names: linking
        # the bare mention would point at a function the reader did not mean.
        nolink = c("catchability", "selectivity")
    ),
    `run-simulation` = list(
        lead   = paste("This guide covers projecting a model forward with",
                       "project(): its arguments, the four ways of specifying",
                       "fishing effort, continuing and comparing runs, and the",
                       "choice of numerical scheme."),
        setup  = c("params <- NS_params", "sim <- NS_sim")
    ),
    `analyse-and-plot` = list(
        # Base generics: a reader clicking these expects base R, not mizer.
        nolink = c("print", "as.data.frame"),
        lead   = paste("This guide gives an overview of the functions",
                       "available in mizer for analysing the results of",
                       "simulations and creating plots. For full documentation",
                       "of each function, follow the links."),
        setup  = c("params <- NS_params", "sim <- NS_sim")
    ),
    `analyse-stability` = list(
        lead   = paste("This guide gives an overview of the",
                       "experimental tools for analysing the dynamic stability",
                       "of a mizer steady state and characterising the limit",
                       "cycles that can replace it. For full documentation of",
                       "each function, follow the links."),
        setup  = "params <- NS_params"
    ),
    `extend-mizer` = list(
        lead   = paste("This guide covers the mechanisms for customising",
                       "mizer's dynamics without editing the package source,",
                       "from external encounter and mortality through to whole",
                       "new ecosystem components, with worked examples of each.",
                       "To package an extension up for others to use, see",
                       "[Creating a mizer extension package](creating-extension-packages.html)."),
        setup  = "params <- NS_params"
    ),
    `use-extension-packages` = list(
        lead   = paste("This guide explains what happens when you load one or",
                       "more mizer extension packages, why the order in which",
                       "you load them can matter, and how to save and share",
                       "models that use extensions. To write an extension",
                       "rather than use one, see the",
                       "[Creating a mizer extension package](creating-extension-packages.html)",
                       "article."),
        # `extensions` is the params slot here, not setExtEncounter()'s page;
        # `library` is base R.
        nolink = c("extensions", "library")
    ),
    # Not a guide, so it overrides all three derived names: the upgrading
    # article is the release-by-release list of changes that break existing
    # code, and the skill an agent loads when a user's script stops working
    # after an upgrade. The skill's symptom index, which is of no use to a
    # human reading by release, is agent-only.
    `upgrade-mizer-code` = list(
        vignette  = "upgrading",
        title     = "Upgrading your mizer code",
        link_text = "Upgrading your mizer code article",
        # Deliberately no `lead` and no `setup`: the skill body opens with its
        # own introduction, and nothing in the article is evaluated.
        nolink    = c("print", "summary", "plot", "as.data.frame")
    )
)

#' Lowercase the first character, leaving the rest alone
#'
#' Turns a heading into something that can follow "the guide to" in a sentence
#' without lowercasing a proper noun later in the phrase.
lcfirst <- function(x) paste0(tolower(substr(x, 1L, 1L)), substr(x, 2L, nchar(x)))

#' The H1 of a skill, which names the topic in prose
skill_h1 <- function(skill, pkg_root = ".") {
    lines <- readLines(file.path(pkg_root, "inst", "skills", skill, "SKILL.md"),
                       warn = FALSE)
    h1 <- grep("^# ", lines, value = TRUE)
    if (!length(h1)) stop("no H1 heading in the ", skill, " skill", call. = FALSE)
    trimws(sub("^# ", "", h1[1]))
}

#' Resolve every topic's three names from its skill
#'
#' Returns a list keyed by skill name with the article's `vignette` basename,
#' its `title`, and the `text` other articles use when linking to it. All three
#' derive from the skill unless the topic overrides them.
guide_index <- function(pkg_root = ".", topics = guide_topics) {
    out <- list()
    for (skill in names(topics)) {
        spec <- topics[[skill]]
        h1 <- skill_h1(skill, pkg_root)
        out[[skill]] <- list(
            vignette = spec$vignette %||% paste0("guide-", skill),
            title    = spec$title %||% paste0("Guide: ", h1),
            text     = spec$link_text %||% paste("guide to", lcfirst(h1))
        )
    }
    out
}


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

#' Read one SKILL.md and return its body: frontmatter and H1 removed.
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

skill_to_guide <- function(skill, spec, index, map, pkg_root = ".") {
    lines <- skill_body(skill, pkg_root)

    # Drop agent-only blocks.
    keep <- rep(TRUE, length(lines))
    inside <- FALSE
    for (i in seq_along(lines)) {
        if (grepl("<!--\\s*agent-only\\s*-->", lines[i])) inside <- TRUE
        if (inside) keep[i] <- FALSE
        if (grepl("<!--\\s*/agent-only\\s*-->", lines[i])) inside <- FALSE
    }
    lines <- lines[keep]

    # Keep article-only blocks but drop their markers: this is the half of the
    # skill that only the article gets, so here it is simply body text.
    # `mizerAgents::setup_mizer_agent()` does the opposite, dropping the block
    # as it installs the skill. The content is knitr source -- a ```{r label}
    # fence inside such a block is evaluated when the article is built, which
    # is the point of keeping it away from an agent, who cannot see the output.
    lines <- lines[!grepl("^\\s*<!--\\s*/?article-only\\s*-->\\s*$", lines)]

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
        for (s in names(index)) {
            target <- index[[s]]
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
    qr_path <- file.path(pkg_root, "inst", "skills", skill,
                         "quick-reference.md")
    if (file.exists(qr_path)) {
        qr <- readLines(qr_path, warn = FALSE)
        qr <- sub("^(\\s*```)r\\s*$", "\\1{r eval=FALSE}", qr)
        out <- c(out, "", "---", "", "## Quick reference", "", qr)
    }

    sources <- paste0("inst/skills/", skill, "/SKILL.md")
    header <- c(
        "---",
        sprintf("title: \"%s\"", index[[skill]]$title),
        "output:",
        "  html_document:",
        "    toc: yes",
        "    toc_float: yes",
        "    fig_width: 5",
        "    fig_height: 4",
        "---",
        "",
        paste0("<!-- Generated from ", sources,
               " by dev_scripts/build_guides.R."),
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
check_output <- function(txt, vignette, index, man_dir, vig_dir) {
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
    is_stray <- vapply(names(index), function(s) {
        grepl(paste0("`", s, "`"), body, fixed = TRUE)
    }, logical(1))
    stray <- names(index)[is_stray]
    if (length(stray)) {
        warning(vignette, ": unconverted skill reference(s): ",
                paste(stray, collapse = ", "),
                ". Phrase each as \"`<skill-name>` skill\" in the source.",
                call. = FALSE)
    }
    length(dead) + length(stray)
}

build_guides <- function(pkg_root = ".",
                         skills = names(guide_topics)) {
    man_dir <- file.path(pkg_root, "man")
    vig_dir <- file.path(pkg_root, "vignettes")
    map <- rd_alias_map(man_dir)
    index <- guide_index(pkg_root)
    bad <- 0L
    for (skill in skills) {
        spec <- guide_topics[[skill]]
        vig <- index[[skill]]$vignette
        txt <- skill_to_guide(skill, spec, index, map, pkg_root)
        bad <- bad + check_output(txt, vig, index, man_dir, vig_dir)
        dest <- file.path(pkg_root, "vignettes", paste0(vig, ".Rmd"))
        writeLines(txt, dest)
        message("Wrote ", dest, " from ", skill)
    }
    if (bad == 0L) {
        message("All links resolve, all skill references converted.")
    }
    invisible(TRUE)
}
