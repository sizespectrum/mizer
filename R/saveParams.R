#' Save and restore mizer objects
#'
#' `saveParams()` saves a MizerParams object to a file. This can then be
#' restored with `readParams()`. `saveSim()` and `readSim()` provide the same
#' lifecycle for MizerSim objects.
#'
#' While these functions ultimately use [saveRDS()] and [readRDS()], they do
#' extra work to make the saved file more robust and more portable, so you
#' should always prefer them over calling `saveRDS()`/`readRDS()` directly on a
#' mizer object.
#'
#' @section What `saveParams()` and `saveSim()` do beyond `saveRDS()`:
#'
#' - They **validate** the object before writing it, so a corrupted or
#'   inconsistent object is caught at save time rather than when you next try
#'   to use it.
#' - They **strip any extension class** and save the object as a plain base
#'   mizer object (recording which extension packages it needs in a slot).
#'   This means the file can be read back even in an R session where the
#'   extension packages that defined those S4 classes are not loaded, and it
#'   protects the file against future changes to those extension classes.
#' - They **check that the required extension packages are installed** and stop
#'   with an informative error if they are not, so you do not save a file that
#'   you would be unable to read back.
#' - They **warn if the model relies on custom functions** (custom rate,
#'   dynamics, selectivity or predation-kernel functions that are not provided
#'   by mizer or a registered extension package). Such functions are not stored
#'   in the file, so to share the model you also need to share an R script or R
#'   Markdown file defining them.
#'
#' Before saving a model you may want to set its metadata with [setMetadata()].
#'
#' @section What `readParams()` and `readSim()` do beyond `readRDS()`:
#'
#' - They **upgrade** an object saved by an older version of mizer to the
#'   current structure (see [upgradeParams()]), so that models saved years ago
#'   still load correctly.
#' - They **re-register the extension packages** that the model needs and,
#'   optionally, install any that are missing (see `install_extensions`),
#'   before restoring the object's extension class.
#' - They **coerce the object back to its extension class** and revalidate it,
#'   reversing the class-stripping done at save time so you get back an object
#'   of the same class you saved.
#'
#' @param params A MizerParams object
#' @param file The name of the file or a connection where the object is saved
#'   to or read from.
#' @param install_extensions Logical. Should [readParams()] or [readSim()]
#'   attempt to install missing extension packages before registering the saved
#'   extension chain?
#' @return `saveParams()` and `saveSim()` return NULL invisibly.
#'   `readParams()` returns a MizerParams object. `readSim()` returns a MizerSim
#'   object.
#' @seealso "Using mizer extension packages":
#'   \code{vignette("using-extension-packages", package = "mizer")}
#' @export
#' @examples
#' # Save params to a temporary file and read them back
#' tmp <- tempfile(fileext = ".rds")
#' saveParams(NS_params, file = tmp)
#' params <- readParams(tmp)
#' identical(params, NS_params)
#'
#' # Save and read back a simulation
#' tmp2 <- tempfile(fileext = ".rds")
#' saveSim(NS_sim, file = tmp2)
#' sim <- readSim(tmp2)
#' identical(sim, NS_sim)
saveParams <- function(params, file) {
    UseMethod("saveParams")
}
#' @export
saveParams.MizerParams <- function(params, file) {
    params <- validateParamsForSaving(params)
    checkRequiredExtensionPackages(params)
    checkCustomFunctions(params)

    saveRDS(baseMizerClass(params), file = file)
}

#' @rdname saveParams
#' @export
readParams <- function(file, install_extensions = FALSE) {
    params <- readRDS(file)
    if (needs_upgrading(params)) {
        params <- suppressWarnings(upgradeParams(params))
    }

    if (length(params@extensions) > 0) {
        registerExtensions(extensionRequirements(params@extensions),
                           install = install_extensions)
    }

    params <- coerceToExtensionClass(params)
    params <- validParams(params)

    # # Check for missing packages
    # packages <- names(params@extensions)
    # missing <- !sapply(packages, require, quietly = TRUE)
    # if (any(missing)) {
    #     warning("Some required extension packages are not installed: ",
    #             paste(missing, collapse = ", "),
    #             ". Shall I install them now? Enter 1 for Yes, ",
    #             " 0 for No.")
    #     ans <- as.integer(readline())
    #     if (ans != 1) return(FALSE)
    #     sapply(packages[missing], remotes::install_github)
    # }
    #
    # # Check for missing functions
    # funs <- c(params@rates_funcs,
    #           params@resource_dynamics,
    #           params@other_dynamics,
    #           params@other_encounter,
    #           params@other_mort,
    #           unique(params@gear_params$sel_func),
    #           paste0(unique(params@species_params$pred_kernel_type),
    #                  "_pred_kernel"))
    # missing <- !sapply(funs, exists, mode = "function")
    # if (any(missing)) {
    #     warning("This model is using the functions ",
    #             paste(funs[missing], collapse = ", "),
    #             ". You need an R script or R Markdown file ",
    #             "defining these functions.")
    # }

    params
}

#' @rdname saveParams
#' @param sim A MizerSim object
#' @export
saveSim <- function(sim, file) {
    UseMethod("saveSim")
}

#' @export
saveSim.MizerSim <- function(sim, file) {
    sim <- validateSimForSaving(sim)
    checkRequiredExtensionPackages(sim@params)
    checkCustomFunctions(sim@params)

    saveRDS(baseMizerClass(sim), file = file)
}

#' @rdname saveParams
#' @export
readSim <- function(file, install_extensions = FALSE) {
    sim <- readRDS(file)
    if (needs_upgrading(sim)) {
        sim <- suppressWarnings(upgradeSim(sim))
    }

    if (length(sim@params@extensions) > 0) {
        registerExtensions(extensionRequirements(sim@params@extensions),
                           install = install_extensions)
    }

    sim@params <- coerceToExtensionClass(sim@params)
    sim <- coerceToExtensionClass(sim)
    sim <- validSim(sim)

    sim
}

validateParamsForSaving <- function(params) {
    tryCatch(
        validParams(params),
        error = function(e) {
            checkRequiredExtensionPackages(params)
            stop(e)
        }
    )
}

validateSimForSaving <- function(sim) {
    tryCatch(
        validSim(sim),
        error = function(e) {
            checkRequiredExtensionPackages(sim@params)
            stop(e)
        }
    )
}

checkRequiredExtensionPackages <- function(params) {
    packages <- requiredExtensionPackages(params)
    missing <- !sapply(packages, requireNamespace, quietly = TRUE)
    if (any(missing)) {
        stop("Some required extension packages are not installed: ",
             paste(packages[missing], collapse = ", "))
    }
}

checkCustomFunctions <- function(params) {
    packages <- requiredExtensionPackages(params)
    kernel_fns <- paste0(unique(params@species_params$pred_kernel_type),
                         "_pred_kernel")
    funs <- c(params@rates_funcs,
              params@resource_dynamics,
              params@other_dynamics,
              params@other_encounter,
              params@other_mort,
              unique(params@gear_params$sel_func),
              kernel_fns)
    custom <- sapply(funs, is_custom, packages = packages)
    if (any(custom)) {
        warning("Your model is using the functions ",
                paste(funs[custom], collapse = ", "),
                ". To share your model you need to share not only the ",
                "params object but also an R script or R Markdown file ",
                "defining these functions.")
    }
}

requiredExtensionPackages <- function(params) {
    reqs <- extensionRequirements(params@extensions)
    c("mizer", names(reqs)[!is.na(reqs)])
}

is_custom <- function(name, packages) {
    !any(sapply(packages,
                function(x, fun) {
                    exists(fun, where = paste0("package:", x),
                           mode = "function")
                },
                fun = name)
         )
}
