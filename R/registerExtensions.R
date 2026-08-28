# Extension tracking and coercion helpers
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later

#' Coerce a mizer object to its extension class
#'
#' Coerces a `MizerParams` or `MizerSim` object to the S3 class hierarchy
#' corresponding to the object's own extension list. For `MizerSim`, the
#' extension names are read from `sim$params$extensions` (or `sim@params@extensions`).
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @param extensions Optional extension list or vector. Defaults to the extensions stored in
#'   `object`, or in `object$params` for `MizerSim`.
#'
#' @return The same object coerced to the appropriate S3 class vector, or to the
#'   base class for an empty extension list.
#' @seealso "Creating a mizer extension package":
#'   \href{https://sizespectrum.org/mizer/articles/guide-create-extension-package.html}{Creating a mizer extension package}
#' @export
#' @family extension tools
coerceToExtensionClass <- function(object, extensions = objectExtensions(object)) {
    if (inherits(object, "MizerParams")) {
        family <- "params"
    } else if (inherits(object, "MizerSim")) {
        family <- "sim"
    } else {
        stop("Can only coerce MizerParams or MizerSim objects.")
    }

    extensions <- validateExtensionsVector(extensions)
    ext_names <- names(extensions)

    if (family == "params") {
        class(object) <- unique(c(ext_names, "MizerParams"))
    } else {
        sim_classes <- if (length(ext_names) > 0) simExtensionClass(ext_names) else character()
        class(object) <- unique(c(sim_classes, "MizerSim"))
    }

    object
}

#' Get the extensions stored in a mizer object
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @return A named character vector of extensions, or an empty character vector
#'   if the object carries no extensions.
#' @keywords internal
objectExtensions <- function(object) {
    if (inherits(object, "MizerParams")) {
        return(extensionRequirements(object$extensions))
    }
    if (inherits(object, "MizerSim")) {
        return(extensionRequirements(object$params$extensions))
    }
    stop("Can only get extensions for MizerParams or MizerSim objects.")
}

#' Extract the requirement view of an extension chain
#'
#' The `$extensions` element may be stored either as a named character vector of
#' requirement strings (the legacy/unversioned form) or as a named list whose
#' entries are length-2 character vectors `c(requirement = ..., version = ...)`.
#' This helper returns the requirement strings as a plain named character
#' vector.
#'
#' @param ext The contents of an `$extensions` element (character vector or list).
#' @return A named character vector of requirement strings.
#' @keywords internal
extensionRequirements <- function(ext) {
    if (is.null(ext) || length(ext) == 0) {
        return(character())
    }
    if (is.list(ext)) {
        req <- vapply(ext, function(e) {
            if (!is.null(names(e)) && "requirement" %in% names(e)) {
                as.character(e[["requirement"]])
            } else {
                as.character(e[[1]])
            }
        }, character(1))
        names(req) <- names(ext)
        return(req)
    }
    ext
}

#' Extract the version stamps of an extension chain
#'
#' Returns the version of the extension package that last upgraded the object
#' for each extension, or `NA_character_` where no stamp is recorded (including
#' the legacy character-vector form, which carries no versions).
#'
#' @param ext The contents of an `$extensions` element (character vector or list).
#' @return A named character vector of version strings (`NA` where unknown).
#' @keywords internal
extensionVersions <- function(ext) {
    if (is.null(ext) || length(ext) == 0) {
        return(character())
    }
    if (is.list(ext)) {
        ver <- vapply(ext, function(e) {
            if (!is.null(names(e)) && "version" %in% names(e)) {
                as.character(e[["version"]])
            } else {
                NA_character_
            }
        }, character(1))
        names(ver) <- names(ext)
        return(ver)
    }
    setNames(rep(NA_character_, length(ext)), names(ext))
}

#' Build a versioned extension list from requirements and versions
#'
#' @param requirements A named character vector of requirement strings.
#' @param versions A named character vector of version stamps. Names not present
#'   default to `NA_character_`.
#' @return A named list whose entries are
#'   `c(requirement = ..., version = ...)`, or an empty character vector when
#'   `requirements` is empty.
#' @keywords internal
makeExtensions <- function(requirements, versions = character()) {
    if (length(requirements) == 0) {
        return(character())
    }
    out <- lapply(names(requirements), function(name) {
        version <- if (name %in% names(versions)) {
            as.character(versions[[name]])
        } else {
            NA_character_
        }
        c(requirement = unname(as.character(requirements[[name]])),
          version = version)
    })
    names(out) <- names(requirements)
    out
}

#' Record an extension and its version stamp on a mizer object
#'
#' Writes an entry for `name` into the object's `$extensions` element, converting
#' it to the versioned list form. Existing entries (and their version stamps)
#' are preserved, keeping their position in the chain. A genuinely new entry is
#' prepended to the front of the list.
#'
#' @param params A `MizerParams` object.
#' @param name The extension identifier (package/class name).
#' @param version Optional version string to stamp. If `NULL` (default) the
#'   existing stamp is preserved.
#' @param requirement Optional requirement string (e.g. `"sizespectrum/mizerMR"`
#'   or `"1.0.0"`).
#' @return The `params` object with the updated `$extensions` element.
#' @seealso "Creating a mizer extension package":
#'   \href{https://sizespectrum.org/mizer/articles/guide-create-extension-package.html}{Creating a mizer extension package}
#' @export
#' @family extension tools
recordExtension <- function(params, name, version = NULL, requirement = NULL) {
    if (isS4(params)) {
        params <- upgrade_s4_to_s3(params)
    }
    assert_that(inherits(params, "MizerParams"), is.string(name))
    ext <- params$extensions
    reqs <- extensionRequirements(ext)
    present <- name %in% names(reqs)

    if (present) {
        req <- if (!is.null(requirement)) as.character(requirement) else unname(reqs[[name]])
    } else {
        req <- if (!is.null(requirement)) as.character(requirement) else NA_character_
    }

    if (is.null(version)) {
        if (present && is.null(requirement)) return(params)
        if (is.list(ext)) {
            entry <- setNames(list(c(requirement = req, version = NA_character_)),
                              name)
            ext <- if (present) { ext[[name]] <- entry[[1]]; ext } else c(entry, ext)
        } else {
            entry <- setNames(req, name)
            ext <- if (present) { ext[[name]] <- req; ext } else c(entry, ext)
        }
        params$extensions <- ext
        return(params)
    }

    if (!is.list(ext)) {
        ext <- makeExtensions(reqs, extensionVersions(ext))
        if (!is.list(ext)) ext <- setNames(list(), character())
    }
    entry <- c(requirement = req, version = as.character(version))
    if (present) {
        ext[[name]] <- entry
    } else {
        ext <- c(setNames(list(entry), name), ext)
    }
    params$extensions <- ext
    params
}

#' Get the recorded version stamp for one extension on an object
#'
#' @param params A `MizerParams` object.
#' @param name The extension identifier.
#' @return The recorded version string, or `NA_character_` if none.
#' @keywords internal
extensionVersion <- function(params, name) {
    vers <- extensionVersions(params$extensions)
    if (name %in% names(vers)) vers[[name]] else NA_character_
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
    if (is.list(extensions)) {
        extensions <- extensionRequirements(extensions)
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
            "Extension names must be syntactically valid class names: ",
            paste(names(extensions)[!valid_names], collapse = ", ")
        )
    }
    extensions
}

#' Derive the MizerSim extension class name for a given extension
#'
#' @param extension Character string — the extension (params) class name.
#' @return A character string formed by appending `"Sim"` to `extension`.
#' @keywords internal
simExtensionClass <- function(extension) {
    paste0(extension, "Sim")
}

#' Load (and optionally install) namespaces for all non-NA extensions
#'
#' @param extensions Named character vector of extensions.
#' @param install Logical. If `TRUE`, install or upgrade missing/outdated
#'   packages via [pak::pkg_install()].
#' @return Invisibly `TRUE`.
#' @keywords internal
ensureExtensionNamespaces <- function(extensions, install = FALSE) {
    extensions <- validateExtensionsVector(extensions)
    if (length(extensions) == 0) {
        return(invisible(TRUE))
    }

    for (extension in names(extensions)) {
        requirement <- unname(extensions[[extension]])

        if (is.na(requirement)) {
            next
        }

        if (!requireNamespace(extension, quietly = TRUE)) {
            if (isTRUE(install)) {
                pkg_spec <- if (isVersionRequirement(requirement)) {
                    extension
                } else {
                    requirement
                }
                pak::pkg_install(pkg_spec)
            }
            if (!requireNamespace(extension, quietly = TRUE)) {
                stop("Required extension package `", extension, "` is not installed.")
            }
        }

        if (isVersionRequirement(requirement) &&
            utils::packageVersion(extension) < package_version(requirement)) {
            if (isTRUE(install)) {
                pak::pkg_install(paste0(extension, "@>=", requirement))
            } else {
                stop(
                    "Extension package `", extension, "` must be at least ",
                    "version ", requirement, ". Use `install = TRUE` to upgrade."
                )
            }
        }

        if (!isNamespaceLoaded(extension)) {
            loadNamespace(extension)
        }
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

#' Strip extension classes from a mizer object
#'
#' Coerces a `MizerParams` or `MizerSim` object back to its plain base class,
#' removing any S3 extension classes. For `MizerSim`, also strips the
#' extension class from the embedded `params` element.
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @return The same object coerced to `MizerParams` or `MizerSim`.
#' @keywords internal
baseMizerClass <- function(object) {
    if (inherits(object, "MizerParams")) {
        class(object) <- "MizerParams"
        object
    } else if (inherits(object, "MizerSim")) {
        class(object) <- "MizerSim"
        if (inherits(object$params, "MizerParams")) {
            class(object$params) <- "MizerParams"
        }
        object
    } else {
        stop("Can only strip extension classes from MizerParams or MizerSim objects.")
    }
}
