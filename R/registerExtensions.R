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

#' Get the extension chain stored in a mizer object
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @return A named character vector of extensions, or an empty character vector
#'   if the object carries no extensions.
#' @keywords internal
objectExtensions <- function(object) {
    if (is(object, "MizerParams")) {
        return(object@extensions)
    }
    if (is(object, "MizerSim")) {
        return(object@params@extensions)
    }
    stop("Can only get extensions for MizerParams or MizerSim objects.")
}

#' Assert that an object's extension chain is compatible with the session
#'
#' Stops with an informative error if the object's extension chain is not a
#' suffix of the session's registered maximal chain, or (when `check_class` is
#' `TRUE`) if the object does not inherit from the expected S4 marker class.
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @param extensions Named character vector giving the object's extension chain.
#'   Defaults to [objectExtensions()] applied to `object`.
#' @param check_class Logical. If `TRUE` (default), also verify that `object`
#'   inherits from the expected S4 marker class.
#' @return Invisibly `TRUE`. Called for its side-effect of stopping on
#'   incompatibility.
#' @keywords internal
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

#' Validate and normalise an extensions named character vector
#'
#' Checks that `extensions` is a named character vector with unique,
#' syntactically valid names, and normalises `NULL` to `character()`.
#'
#' @param extensions A named character vector, or `NULL`.
#' @return A validated named character vector (possibly length-zero).
#' @keywords internal
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

#' Compare two extension chains
#'
#' @param old Named character vector for the previously registered chain.
#' @param new Named character vector for the proposed chain.
#' @return One of `"identical"`, `"new_is_suffix"`, `"old_is_suffix"`, or
#'   `"incompatible"`.
#' @keywords internal
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

#' Test whether one extension chain is a suffix of another
#'
#' An empty `candidate` is always a suffix. Order and values must match
#' exactly for the overlapping tail.
#'
#' @param candidate Named character vector to test.
#' @param chain Named character vector that may contain `candidate` as a tail.
#' @return `TRUE` if `candidate` is a suffix of `chain`, `FALSE` otherwise.
#' @keywords internal
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

#' Define S4 marker classes for a set of dispatch extensions
#'
#' Creates a linear inheritance chain of S4 classes: the outermost extension
#' extends the next, which extends the next, down to the base `MizerParams` /
#' `MizerSim` class. Existing classes are checked for compatibility instead of
#' being redefined.
#'
#' @param extensions Named character vector of extensions (full chain or
#'   dispatch subset). Non-dispatch entries are silently ignored.
#' @return Invisibly, the named character vector of dispatch extensions.
#' @keywords internal
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

#' Filter an extension vector to those that participate in S3/S4 dispatch
#'
#' An extension participates in dispatch if its requirement is `NA_character_`
#' (in-development) or if an S4 class with its name already exists.
#'
#' @param extensions Named character vector of extensions.
#' @return A named character vector containing only the dispatch extensions,
#'   preserving order.
#' @keywords internal
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

#' Test whether a mizer object uses extension S4 dispatch
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @return `TRUE` if the object's primary class is not the plain base class.
#' @keywords internal
usesExtensionDispatch <- function(object) {
    if (is(object, "MizerParams")) {
        return(!identical(class(object)[[1]], "MizerParams"))
    }
    if (is(object, "MizerSim")) {
        return(!identical(class(object)[[1]], "MizerSim"))
    }
    stop("Can only check dispatch for MizerParams or MizerSim objects.")
}

#' Define an S4 class or verify it extends the expected parent
#'
#' If `class` does not yet exist, defines it as a virtual-free S4 class that
#' contains `parent`, registered in `.GlobalEnv`. If `class` already exists,
#' stops with an error unless it already extends `parent`.
#'
#' @param class Character string — the S4 class name to define or check.
#' @param parent Character string — the required parent class.
#' @return Invisibly, `class`.
#' @keywords internal
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

#' Derive the MizerSim marker class name for a given extension
#'
#' @param extension Character string — the extension (params) class name.
#' @return A character string formed by appending `"Sim"` to `extension`.
#' @keywords internal
simExtensionClass <- function(extension) {
    paste0(extension, "Sim")
}

#' Load (and optionally install) namespaces for all non-NA extensions
#'
#' For each extension whose requirement is not `NA_character_`, checks that the
#' package is installed, installs it if `install = TRUE` and it is missing,
#' verifies the minimum version if the requirement is a version string, then
#' calls [loadNamespace()].
#'
#' @param extensions Named character vector of extensions.
#' @param install Logical. If `TRUE`, attempt to install missing packages via
#'   [utils::install.packages()].
#' @return Invisibly `TRUE`.
#' @keywords internal
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

#' Test whether a requirement string is a dotted version number
#'
#' @param requirement Character string.
#' @return `TRUE` if `requirement` matches `"X.Y.Z..."` (digits and dots only).
#' @keywords internal
isVersionRequirement <- function(requirement) {
    grepl("^[0-9]+(\\.[0-9]+)*$", requirement)
}

#' Format an extension chain as a human-readable string
#'
#' @param extensions Named character vector of extensions.
#' @return A character string such as `"mizerExtB -> mizerExtA"`, or
#'   `"<empty>"` for a zero-length chain.
#' @keywords internal
formatExtensionChain <- function(extensions) {
    if (length(extensions) == 0) {
        return("<empty>")
    }
    paste(names(extensions), collapse = " -> ")
}

#' Strip extension classes from a mizer object
#'
#' Coerces a `MizerParams` or `MizerSim` object back to its plain base class,
#' removing any S4 extension marker classes. For `MizerSim`, also strips the
#' extension class from the embedded `params` slot.
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @return The same object coerced to `MizerParams` or `MizerSim`.
#' @keywords internal
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

#' Reset the registered extension chain for the current R session
#'
#' Clears the session's extension registry. Primarily used in tests to restore
#' a clean state between test cases.
#'
#' @return Invisibly, an empty character vector.
#' @keywords internal
resetMizerSession <- function() {
    .mizerSession$extensions <- character()
    invisible(character())
}
