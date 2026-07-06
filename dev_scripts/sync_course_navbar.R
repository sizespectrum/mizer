#!/usr/bin/env Rscript

# Synchronise the mizerCourse Quarto top navbar with the mizer pkgdown navbar.
#
# Usage:
#   Rscript dev_scripts/sync_course_navbar.R
#   Rscript dev_scripts/sync_course_navbar.R --course ../mizerCourse
#
# The script intentionally replaces only the `website.navbar` block in the
# course _quarto-website.yml, leaving the sidebar and the rest of the Quarto
# configuration untouched.

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) {
        y
    } else {
        x
    }
}

value_after <- function(flag, default = NULL) {
    pos <- match(flag, args)
    if (is.na(pos) || pos == length(args)) {
        return(default)
    }
    args[[pos + 1]]
}

get_script_path <- function() {
    file_arg <- grep(
        "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
    )
    if (length(file_arg) > 0) {
        return(sub("^--file=", "", file_arg[1]))
    }
    for (i in rev(seq_len(sys.nframe()))) {
        frame <- sys.frame(i)
        if (exists("ofile", envir = frame, inherits = FALSE)) {
            return(get("ofile", envir = frame))
        }
    }
    "dev_scripts/sync_course_navbar.R"
}

script_path <- get_script_path()

repo_root <- normalizePath(
    file.path(dirname(script_path), ".."),
    mustWork = TRUE
)

pkgdown_yml <- normalizePath(
    value_after("--pkgdown", file.path(repo_root, "pkgdown", "_pkgdown.yml")),
    mustWork = TRUE
)
course_dir <- normalizePath(
    value_after("--course", file.path(repo_root, "..", "mizerCourse")),
    mustWork = TRUE
)
course_yml <- normalizePath(
    value_after("--quarto", file.path(course_dir, "_quarto-website.yml")),
    mustWork = TRUE
)

need_package <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Please install the '", package, "' package.", call. = FALSE)
    }
}

need_package("yaml")

read_config <- function(path) {
    yaml::read_yaml(path, eval.expr = FALSE)
}

yaml_null <- function() {
    structure(list(), class = "yaml_null")
}

is_yaml_null <- function(x) {
    is.null(x) || inherits(x, "yaml_null")
}

site_base_url <- function(config) {
    sub("/+$", "", config$url %||% "https://sizespectrum.org/mizer")
}

site_href <- function(path, base_url) {
    if (grepl("^https?://", path)) {
        return(path)
    }
    paste0(base_url, "/", sub("^/+", "", path))
}

article_href <- function(slug, base_url) {
    site_href(file.path("articles", paste0(slug, ".html")), base_url)
}

pretty_label <- function(slug) {
    words <- gsub("[-_]+", " ", slug)
    tools::toTitleCase(words)
}

read_source_title <- function(slug, repo_root) {
    candidates <- file.path(repo_root, "vignettes", paste0(slug, c(".Rmd", ".qmd")))
    path <- candidates[file.exists(candidates)][1]
    if (is.na(path)) {
        return(NULL)
    }

    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    if (!identical(trimws(lines[[1]] %||% ""), "---")) {
        return(NULL)
    }

    end <- which(trimws(lines[-1]) == "---")[1]
    if (is.na(end)) {
        return(NULL)
    }
    front_matter <- paste(lines[seq.int(2, end)], collapse = "\n")
    title <- tryCatch(
        yaml::yaml.load(front_matter, eval.expr = FALSE)$title,
        error = function(e) NULL
    )
    if (is.null(title)) {
        return(NULL)
    }
    as.character(title)
}

read_article_title <- function(slug, repo_root) {
    source_title <- read_source_title(slug, repo_root)
    if (!is.null(source_title) && nzchar(source_title)) {
        return(source_title)
    }

    path <- file.path(repo_root, "docs", "articles", paste0(slug, ".html"))
    if (!file.exists(path)) {
        return(pretty_label(slug))
    }

    html <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    title <- sub("(?is).*?<h1[^>]*>(.*?)</h1>.*", "\\1", html, perl = TRUE)
    if (identical(title, html)) {
        title <- sub("(?is).*<title[^>]*>(.*?)</title>.*", "\\1", html, perl = TRUE)
    }
    if (identical(title, html)) {
        return(pretty_label(slug))
    }

    title <- gsub("<[^>]+>", "", title)
    title <- gsub("\\s+", " ", title)
    title <- trimws(title)
    title <- sub(" \\| mizer$", "", title)
    if (nzchar(title)) title else pretty_label(slug)
}

normalise_contents <- function(contents) {
    if (is.null(contents)) {
        return(character())
    }
    contents <- unlist(contents, use.names = FALSE)
    contents <- contents[!grepl("^(has_concept|starts_with)\\(", contents)]
    contents <- contents[!startsWith(contents, "-")]
    contents
}

build_articles_menu <- function(config, base_url, repo_root, intro_slug) {
    menu <- list()

    for (section in config$articles %||% list()) {
        if (!"navbar" %in% names(section)) {
            next
        }
        contents <- normalise_contents(section$contents)
        contents <- setdiff(contents, intro_slug)
        if (!length(contents)) {
            next
        }

        navbar_label <- section$navbar
        if (!is_yaml_null(navbar_label)) {
            menu[[length(menu) + 1]] <- list(text = as.character(navbar_label))
        }

        for (slug in contents) {
            menu[[length(menu) + 1]] <- list(
                text = read_article_title(slug, repo_root),
                href = article_href(slug, base_url)
            )
        }
    }

    menu[[length(menu) + 1]] <- list(
        text = "More articles...",
        href = site_href("articles/index.html", base_url)
    )
    menu
}

intro_slug <- function(config) {
    articles <- unlist(lapply(config$articles %||% list(), function(section) {
        normalise_contents(section$contents)
    }), use.names = FALSE)

    rproj <- list.files(
        repo_root, pattern = "\\.Rproj$", full.names = FALSE
    )
    package <- sub("\\.Rproj$", "", basename(rproj[1] %||% "mizer"))

    package_vignette <- file.path(repo_root, "vignettes", paste0(package, c(".Rmd", ".qmd")))

    if (any(file.exists(package_vignette)) || package %in% articles) {
        package
    } else {
        articles[[1]] %||% package
    }
}

component_to_quarto <- function(name, component, config, base_url, repo_root, intro) {
    builtin <- switch(name,
        home = NULL,
        intro = list(
            text = "Get started",
            href = article_href(intro, base_url)
        ),
        articles = list(
            text = "Articles",
            menu = build_articles_menu(config, base_url, repo_root, intro)
        ),
        reference = list(
            text = "Reference",
            href = site_href("reference/index.html", base_url)
        ),
        news = list(
            text = "Changelog",
            href = site_href("news/index.html", base_url)
        ),
        search = NULL,
        lightswitch = NULL,
        github = list(
            icon = "github",
            href = "https://github.com/sizespectrum/mizer/",
            `aria-label` = "GitHub"
        ),
        NULL
    )
    if (!is.null(builtin) || name %in% c("home", "search", "lightswitch")) {
        return(builtin)
    }

    if (is_yaml_null(component)) {
        return(NULL)
    }

    if (is.null(component$text) && is.null(component$href) && is.null(component$icon)) {
        stop("Do not know how to translate pkgdown navbar component '", name, "'.",
             call. = FALSE)
    }

    item <- list()
    if (!is.null(component$icon)) {
        icon <- sub("^fa-", "", component$icon)
        icon <- sub("^fab fa-", "", icon)
        icon <- sub(" fa-lg$", "", icon)
        item$icon <- icon
    }
    if (!is.null(component$text)) {
        item$text <- component$text
    }
    if (!is.null(component$href)) {
        item$href <- site_href(component$href, base_url)
    }
    if (!is.null(component$`aria-label`)) {
        item$`aria-label` <- component$`aria-label`
    }
    item
}

build_nav_side <- function(names, components, config, base_url, repo_root, intro) {
    items <- list()
    for (name in names %||% character()) {
        component <- components[[name]]
        item <- component_to_quarto(name, component, config, base_url, repo_root, intro)
        if (!is.null(item)) {
            items[[length(items) + 1]] <- item
        }
    }
    items
}

build_navbar <- function(config, repo_root) {
    base_url <- site_base_url(config)
    navbar <- config$navbar
    structure <- navbar$structure
    components <- navbar$components %||% list()
    intro <- intro_slug(config)

    list(
        title = "mizer",
        `logo-href` = base_url,
        search = "search" %in% unlist(structure, use.names = FALSE),
        left = build_nav_side(
            structure$left, components, config, base_url, repo_root, intro
        ),
        right = build_nav_side(
            structure$right, components, config, base_url, repo_root, intro
        )
    )
}

as_yaml_block <- function(navbar) {
    lines <- yaml::as.yaml(
        list(navbar = navbar),
        indent = 2,
        indent.mapping.sequence = TRUE,
        line.sep = "\n",
        handlers = list(yaml_null = function(x) "null")
    )
    lines <- strsplit(lines, "\n", fixed = TRUE)[[1]]
    lines <- sub(": yes$", ": true", lines)
    lines <- sub(": no$", ": false", lines)
    paste0("  ", lines)
}

replace_navbar_block <- function(lines, block) {
    website_line <- grep("^website:\\s*$", lines)
    if (length(website_line) != 1) {
        stop("Could not find a unique top-level website: block.", call. = FALSE)
    }

    start <- grep("^  navbar:\\s*$", lines)
    start <- start[start > website_line]
    if (!length(start)) {
        start <- grep("^  page-navigation:", lines)
        start <- start[start > website_line]
        if (!length(start)) {
            stop("Could not find where to insert website.navbar.", call. = FALSE)
        }
        start <- start[[1]]
        return(append(lines, block, after = start - 1))
    }
    start <- start[[1]]

    end <- start + 1
    while (end <= length(lines) &&
           (grepl("^    ", lines[[end]]) || grepl("^\\s*$", lines[[end]]))) {
        end <- end + 1
    }

    after_block <- if (end <= length(lines)) lines[end:length(lines)] else character()
    c(lines[seq_len(start - 1)], block, after_block)
}

config <- read_config(pkgdown_yml)
navbar <- build_navbar(config, repo_root)
navbar_block <- as_yaml_block(navbar)

course_lines <- readLines(course_yml, warn = FALSE, encoding = "UTF-8")
updated_lines <- replace_navbar_block(course_lines, navbar_block)

if (identical(course_lines, updated_lines)) {
    message("Course navbar is already up to date: ", course_yml)
} else {
    writeLines(updated_lines, course_yml, useBytes = TRUE)
    message("Updated course navbar: ", course_yml)
}
