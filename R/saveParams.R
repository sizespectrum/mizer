#' Save and restore mizer objects
#'
#' `saveParams()` saves a MizerParams object to a file. This can then be
#' restored with `readParams()`. `saveSim()` and `readSim()` provide the same
#' lifecycle for MizerSim objects.
#'
#' Issues a warning if the model you are saving relies on some custom functions.
#' Before saving a model you may want to set its metadata with
#' [setMetadata()].
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
        registerExtensions(params@extensions, install = install_extensions)
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
        registerExtensions(sim@params@extensions, install = install_extensions)
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
    c("mizer", names(params@extensions)[!is.na(params@extensions)])
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
