# Extension chain registration and coercion helpers
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later

.mizerSession <- new.env(parent = emptyenv())
.mizerSession$extensions <- character()

#' Register a single mizer extension for this R session
#'
#' Prepends one extension to the front of the active extension chain, giving it
#' the highest dispatch priority. Designed to be called from a package's
#' `.onLoad` hook so that the chain grows naturally in load order: the last
#' package loaded ends up outermost.
#'
#' The call is idempotent: if the extension is already registered at any
#' position in the chain, the function leaves the chain unchanged, including the
#' requirement recorded at first registration.
#'
#' @param name A syntactically valid R name identifying the extension (e.g.
#'   `"mizerExtA"`). This name is used as the S3 marker class name.
#' @param requirement A version string, installation specification, or
#'   `NA_character_` (the default). `NA_character_` marks an in-development
#'   extension. A version string such as `"1.2.0"` records the minimum required
#'   package version.
#' @param install Logical. If `TRUE`, attempt to install a missing extension
#'   package.
#'
#' @return The updated extension chain, invisibly.
#' @seealso [registerExtensions()] for registering an explicit full chain.
#'   The [guide to using mizer extension
#'   packages](https://sizespectrum.org/mizer/articles/guide-use-extension-packages.html).
#'   "Creating a mizer extension package":
#'   \href{https://sizespectrum.org/mizer/articles/guide-create-extension-package.html}{Creating a mizer extension package}
#' @export
#' @family extension tools
registerExtension <- function(name, requirement = NA_character_, install = FALSE) {
    extension <- validateExtensionsVector(setNames(requirement, name))

    current <- getRegisteredExtensions()

    if (name %in% names(current)) {
        return(invisible(current))
    }

    new_chain <- c(extension, current)
    ensureExtensionNamespaces(new_chain, install = install)
    .mizerSession$extensions <- new_chain

    invisible(new_chain)
}

#' Register mizer extensions for this R session
#'
#' Registers an explicit full extension chain for the current R session. The
#' order of `extensions` is the S3 dispatch order, from outermost to innermost
#' extension. For example `c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")`
#' dispatches to `mizerExtB` methods first, then `mizerExtA` methods, then base
#' mizer methods.
#'
#' A session can handle objects whose extension chain is a suffix of the
#' registered maximal chain. For example, after registering
#' `c(mizerExtB = "1.2.0", mizerExtA = "0.4.1")`, objects using only
#' `c(mizerExtA = "0.4.1")` are also valid.
#'
#' For extension packages that register themselves incrementally from `.onLoad`,
#' use [registerExtension()] instead.
#'
#' @param extensions A named character vector. Names are extension identifiers.
#'   Values are version strings, installation specifications, or
#'   `NA_character_`. Installed extensions participate in S3 dispatch if they
#'   provide S3 methods for their class name. `NA_character_` entries are
#'   treated as in-development dispatch extensions.
#' @param install Logical. If `TRUE`, missing or outdated extension packages are
#'   installed via [pak::pkg_install()]. Version strings install from CRAN;
#'   other requirement strings (e.g. `"user/repo@v1.2.0"`) are passed directly
#'   to pak and may refer to GitHub, local paths, or any other pak-supported
#'   source.
#'
#' @return The active maximal extension chain, invisibly.
#' @seealso [registerExtension()] for the incremental per-package variant.
#'   The [guide to using mizer extension
#'   packages](https://sizespectrum.org/mizer/articles/guide-use-extension-packages.html).
#'   "Creating a mizer extension package":
#'   \href{https://sizespectrum.org/mizer/articles/guide-create-extension-package.html}{Creating a mizer extension package}
#' @export
#' @family extension tools
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
#' @seealso The [guide to using mizer extension
#'   packages](https://sizespectrum.org/mizer/articles/guide-use-extension-packages.html)
#' @export
#' @family extension tools
getRegisteredExtensions <- function() {
    .mizerSession$extensions
}

#' Coerce a mizer object to its registered extension class
#'
#' Coerces a `MizerParams` or `MizerSim` object to the S3 class hierarchy
#' corresponding to the object's own extension chain. For `MizerSim`, the
#' extension chain is read from `sim$params$extensions` (or `sim@params@extensions`).
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @param extensions Optional extension chain. Defaults to the chain stored in
#'   `object`, or in `object$params` for `MizerSim`.
#'
#' @return The same object coerced to the appropriate S3 class vector, or to the
#'   base class for an empty extension chain.
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
    assertExtensionChain(object, extensions = extensions, check_class = FALSE)

    dispatch_extensions <- dispatchExtensions(extensions)

    if (family == "params") {
        class(object) <- c(names(dispatch_extensions), "MizerParams")
    } else {
        sim_classes <- if (length(dispatch_extensions) > 0) paste0(names(dispatch_extensions), "Sim") else character()
        class(object) <- c(sim_classes, "MizerSim")
    }

    object
}

#' Get the extension chain stored in a mizer object
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
#' vector, which is the form used for dispatch and suffix comparison.
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
#' prepended to the front of the chain so that it stays ordered outermost-first,
#' matching [registerExtension()]. The requirement is taken from the existing
#' entry if present, otherwise from the registered extension chain.
#'
#' @param params A `MizerParams` object.
#' @param name The extension identifier.
#' @param version Optional version string to stamp. If `NULL` (default) the
#'   existing stamp is preserved.
#' @return The `params` object with the updated `$extensions` element.
#' @seealso "Creating a mizer extension package":
#'   \href{https://sizespectrum.org/mizer/articles/guide-create-extension-package.html}{Creating a mizer extension package}
#' @export
#' @family extension tools
recordExtension <- function(params, name, version = NULL) {
    if (isS4(params)) {
        params <- upgrade_s4_to_s3(params)
    }
    assert_that(inherits(params, "MizerParams"), is.string(name))
    ext <- params$extensions
    reqs <- extensionRequirements(ext)
    present <- name %in% names(reqs)

    # Requirement: keep the existing one, else take it from the registry.
    if (present) {
        req <- unname(reqs[[name]])
    } else {
        registered <- getRegisteredExtensions()
        req <- if (name %in% names(registered)) {
            unname(registered[[name]])
        } else {
            NA_character_
        }
    }

    if (is.null(version)) {
        if (present) return(params)
        if (is.list(ext)) {
            entry <- setNames(list(c(requirement = req, version = NA_character_)),
                              name)
            ext <- c(entry, ext)
        } else {
            ext <- c(setNames(req, name), ext)
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

#' Assert that an object's extension chain is compatible with the session
#'
#' Stops with an informative error if the object's extension chain is not a
#' suffix of the session's registered maximal chain, or (when `check_class` is
#' `TRUE`) if the object does not inherit from the expected class.
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @param extensions Named character vector giving the object's extension chain.
#'   Defaults to [objectExtensions()] applied to `object`.
#' @param check_class Logical. If `TRUE` (default), also verify that `object`
#'   inherits from the expected class.
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

    if (inherits(object, "MizerParams")) {
        expected_class <- if (length(dispatch_extensions) == 0) {
            "MizerParams"
        } else {
            names(dispatch_extensions)[1]
        }
    } else if (inherits(object, "MizerSim")) {
        expected_class <- if (length(dispatch_extensions) == 0) {
            "MizerSim"
        } else {
            simExtensionClass(names(dispatch_extensions)[1])
        }
    } else {
        stop("Can only check MizerParams or MizerSim objects.")
    }

    if (isTRUE(check_class) && !inherits(object, expected_class)) {
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

#' Does an installed extension register dispatch methods for its own class?
#'
#' @param name The extension identifier (its class / package name).
#' @return `TRUE` if the loaded namespace `name` registers S3 methods for class
#'   `name` or `paste0(name, "Sim")`, otherwise `FALSE`.
#' @keywords internal
providesDispatchMethods <- function(name) {
    if (!isNamespaceLoaded(name)) {
        return(FALSE)
    }
    s3 <- tryCatch(getNamespaceInfo(asNamespace(name), "S3methods"),
                   error = function(e) NULL)
    if (is.null(s3) || nrow(s3) == 0) {
        return(FALSE)
    }
    classes <- s3[, 2]
    name %in% classes || simExtensionClass(name) %in% classes
}

#' Filter an extension vector to those that participate in S3 dispatch
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

        is.na(requirement) || providesDispatchMethods(extension)
    }, logical(1))

    extensions[is_dispatch_extension]
}

#' Test whether a mizer object uses extension S3 dispatch
#'
#' @param object A `MizerParams` or `MizerSim` object.
#' @return `TRUE` if the object's primary class is not the plain base class.
#' @keywords internal
usesExtensionDispatch <- function(object) {
    if (inherits(object, "MizerParams")) {
        return(!identical(class(object)[[1]], "MizerParams"))
    }
    if (inherits(object, "MizerSim")) {
        return(!identical(class(object)[[1]], "MizerSim"))
    }
    stop("Can only check dispatch for MizerParams or MizerSim objects.")
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

#' Clear the registered extension chain
#'
#' Clears the session's extension registry. You can then create a new
#' extension chain with [registerExtensions()].
#'
#' @return Invisibly, an empty character vector.
#' @seealso The [guide to using mizer extension
#'   packages](https://sizespectrum.org/mizer/articles/guide-use-extension-packages.html)
#' @family extension tools
#' @export
clearExtensionChain <- function() {
    .mizerSession$extensions <- character()
    invisible(character())
}
