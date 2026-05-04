# Extension chain registration and coercion helpers
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later

.mizerSession <- new.env(parent = emptyenv())
.mizerSession$extensions <- character()

#' Register mizer extensions for this R session
#'
#' Registers the maximal extension chain that mizer should use in the current
#' R session. The order of `extensions` is the S3 dispatch order, from
#' outermost to innermost extension. For example
#' `c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")` dispatches to
#' `mizerExtB` methods first, then `mizerExtA` methods, then base mizer
#' methods.
#'
#' A session can handle objects whose extension chain is a suffix of the
#' registered maximal chain. For example, after registering
#' `c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")`, objects using only
#' `c(mizerExtA = "0.4.1")` are also valid.
#'
#' @param extensions A named character vector. Names are extension identifiers.
#'   Values are version strings, installation specifications, or
#'   `NA_character_`. Installed extensions only participate in S3 dispatch if
#'   they provide an S4 marker class with the same name. `NA_character_` entries
#'   are treated as in-development dispatch extensions and mizer creates their
#'   marker classes automatically.
#' @param install Logical. If `TRUE`, missing extension packages may be
#'   installed. Installation support is intentionally conservative and currently
#'   only supports CRAN-style package installation by extension name.
#'
#' @return The active maximal extension chain, invisibly.
#' @export
registerExtensions <- function(extensions, install = FALSE) {
    extensions <- validateExtensionsVector(extensions)
    old <- getRegisteredExtensions()
    relation <- compareExtensionChains(old, extensions)

    if (relation %in% c("identical", "new_is_suffix")) {
        ensureExtensionNamespaces(extensions, install = install)
        return(invisible(old))
    }

    if (relation == "old_is_suffix") {
        ensureExtensionNamespaces(extensions, install = install)
        defineExtensionClasses(extensions)
        .mizerSession$extensions <- extensions
        return(invisible(extensions))
    }

    stop(
        "A different extension chain is already active in this session. ",
        "Please restart R before registering this chain."
    )
}

#' Get the registered mizer extension chain
#'
#' @return A named character vector giving the maximal extension chain
#'   registered for this R session.
#' @export
getRegisteredExtensions <- function() {
    .mizerSession$extensions
}

#' Coerce a mizer object to its registered extension class
#'
#' Coerces a `MizerParams` or `MizerSim` object to the S4 marker class
#' corresponding to the object's own extension chain. For `MizerSim`, the
#' extension chain is read from `sim@params@extensions`.
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @param extensions Optional extension chain. Defaults to the chain stored in
#'   `object`, or in `object@params` for `MizerSim`.
#'
#' @return The same object coerced to the appropriate marker class, or to the
#'   base class for an empty extension chain.
#' @export
coerceToExtensionClass <- function(object, extensions = objectExtensions(object)) {
    if (is(object, "MizerParams")) {
        family <- "params"
        base_class <- "MizerParams"
    } else if (is(object, "MizerSim")) {
        family <- "sim"
        base_class <- "MizerSim"
    } else {
        stop("Can only coerce MizerParams or MizerSim objects.")
    }

    extensions <- validateExtensionsVector(extensions)
    assertExtensionChain(object, extensions = extensions, check_class = FALSE)

    dispatch_extensions <- dispatchExtensions(extensions)

    if (length(dispatch_extensions) == 0) {
        return(methods::as(object, base_class))
    }

    target_class <- names(dispatch_extensions)[1]
    if (family == "sim") {
        target_class <- simExtensionClass(target_class)
    }

    methods::as(object, target_class)
}

objectExtensions <- function(object) {
    if (is(object, "MizerParams")) {
        return(object@extensions)
    }
    if (is(object, "MizerSim")) {
        return(object@params@extensions)
    }
    stop("Can only get extensions for MizerParams or MizerSim objects.")
}

assertExtensionChain <- function(object, extensions = objectExtensions(object),
                                 check_class = TRUE) {
    extensions <- validateExtensionsVector(extensions)
    active <- getRegisteredExtensions()

    if (!isSuffixChain(extensions, active)) {
        if (length(extensions) > 0 && length(active) == 0) {
            stop(
                "This object uses mizer extensions but no compatible ",
                "extension chain is registered. Please call ",
                "`registerExtensions()` first, or load the object with ",
                "`readParams()`."
            )
        }
        stop(
            "This object's extension chain is incompatible with the ",
            "extension chain registered in this R session. Please restart R ",
            "before using this object."
        )
    }

    dispatch_extensions <- dispatchExtensions(extensions)

    if (is(object, "MizerParams")) {
        expected_class <- if (length(dispatch_extensions) == 0) {
            "MizerParams"
        } else {
            names(dispatch_extensions)[1]
        }
    } else if (is(object, "MizerSim")) {
        expected_class <- if (length(dispatch_extensions) == 0) {
            "MizerSim"
        } else {
            simExtensionClass(names(dispatch_extensions)[1])
        }
    } else {
        stop("Can only check MizerParams or MizerSim objects.")
    }

    if (isTRUE(check_class) && !is(object, expected_class)) {
        stop(
            "This object has extension chain ",
            formatExtensionChain(extensions),
            " but does not inherit from class `", expected_class, "`."
        )
    }

    invisible(TRUE)
}

validateExtensionsVector <- function(extensions) {
    if (is.null(extensions)) {
        extensions <- character()
    }
    if (!is.character(extensions)) {
        stop("`extensions` must be a named character vector.")
    }
    if (length(extensions) == 0) {
        return(character())
    }
    if (is.null(names(extensions))) {
        stop("`extensions` must be a named character vector.")
    }
    if (anyNA(names(extensions)) || any(names(extensions) == "")) {
        stop("All entries in `extensions` must have non-empty names.")
    }
    if (anyDuplicated(names(extensions))) {
        stop("Extension names must be unique.")
    }
    valid_names <- make.names(names(extensions)) == names(extensions)
    if (!all(valid_names)) {
        stop(
            "Extension names must be syntactically valid S4 class names: ",
            paste(names(extensions)[!valid_names], collapse = ", ")
        )
    }
    extensions
}

compareExtensionChains <- function(old, new) {
    if (identical(old, new)) {
        return("identical")
    }
    if (isSuffixChain(new, old)) {
        return("new_is_suffix")
    }
    if (isSuffixChain(old, new)) {
        return("old_is_suffix")
    }
    "incompatible"
}

isSuffixChain <- function(candidate, chain) {
    candidate <- validateExtensionsVector(candidate)
    chain <- validateExtensionsVector(chain)

    if (length(candidate) == 0) {
        return(TRUE)
    }
    if (length(candidate) > length(chain)) {
        return(FALSE)
    }

    start <- length(chain) - length(candidate) + 1
    identical(candidate, chain[start:length(chain)])
}

defineExtensionClasses <- function(extensions) {
    extensions <- validateExtensionsVector(extensions)
    extensions <- dispatchExtensions(extensions)
    parent_params <- "MizerParams"
    parent_sim <- "MizerSim"

    for (extension in rev(names(extensions))) {
        defineOrCheckClass(extension, parent_params)
        parent_params <- extension

        sim_class <- simExtensionClass(extension)
        defineOrCheckClass(sim_class, parent_sim)
        parent_sim <- sim_class
    }

    invisible(extensions)
}

dispatchExtensions <- function(extensions) {
    extensions <- validateExtensionsVector(extensions)
    if (length(extensions) == 0) {
        return(character())
    }

    is_dispatch_extension <- vapply(seq_along(extensions), function(i) {
        extension <- names(extensions)[[i]]
        requirement <- unname(extensions[[i]])

        is.na(requirement) || methods::isClass(extension)
    }, logical(1))

    extensions[is_dispatch_extension]
}

usesExtensionDispatch <- function(object) {
    if (is(object, "MizerParams")) {
        return(!identical(class(object)[[1]], "MizerParams"))
    }
    if (is(object, "MizerSim")) {
        return(!identical(class(object)[[1]], "MizerSim"))
    }
    stop("Can only check dispatch for MizerParams or MizerSim objects.")
}

defineOrCheckClass <- function(class, parent) {
    if (!methods::isClass(class)) {
        methods::setClass(class, contains = parent, where = .GlobalEnv)
        return(invisible(class))
    }

    if (!methods::extends(class, parent)) {
        stop(
            "Class `", class, "` already exists but does not contain `",
            parent, "`."
        )
    }

    invisible(class)
}

simExtensionClass <- function(extension) {
    paste0(extension, "Sim")
}

ensureExtensionNamespaces <- function(extensions, install = FALSE) {
    extensions <- validateExtensionsVector(extensions)
    if (length(extensions) == 0) {
        return(invisible(TRUE))
    }

    for (extension in names(extensions)) {
        requirement <- unname(extensions[[extension]])

        # NA means "marker class only" for tests or in-development extensions.
        if (is.na(requirement)) {
            next
        }

        if (!requireNamespace(extension, quietly = TRUE)) {
            if (isTRUE(install)) {
                utils::install.packages(extension)
            }
            if (!requireNamespace(extension, quietly = TRUE)) {
                stop("Required extension package `", extension, "` is not installed.")
            }
        }

        if (isVersionRequirement(requirement) &&
            utils::packageVersion(extension) < package_version(requirement)) {
            stop(
                "Extension package `", extension, "` must be at least ",
                "version ", requirement, "."
            )
        }

        loadNamespace(extension)
    }

    invisible(TRUE)
}

isVersionRequirement <- function(requirement) {
    grepl("^[0-9]+(\\.[0-9]+)*$", requirement)
}

formatExtensionChain <- function(extensions) {
    if (length(extensions) == 0) {
        return("<empty>")
    }
    paste(names(extensions), collapse = " -> ")
}

baseMizerClass <- function(object) {
    if (is(object, "MizerParams")) {
        methods::as(object, "MizerParams")
    } else if (is(object, "MizerSim")) {
        sim <- methods::as(object, "MizerSim")
        sim@params <- methods::as(sim@params, "MizerParams")
        sim
    } else {
        stop("Can only strip extension classes from MizerParams or MizerSim objects.")
    }
}

resetMizerSession <- function() {
    .mizerSession$extensions <- character()
    invisible(character())
}
