#' Save a MizerParams object
#' 
#' saveParams()` saves a MizerParams object to a file. This can then be loaded
#' back with `readParams()`.
#' 
#' Issues a warning if the model you are saving relies on some custom functions.
#' 
#' @param params A MizerParams object
#' @param file The name of the file or a connection where the MizerParams object
#'   is saved to or read from.
#' @export
saveParams <- function(params, file) {
    params <- validParams(params)
    
    kernel_fns <- paste0(unique(params@species_params$pred_kernel_type),
                         "_pred_kernel")
    funs <- c(params@rates_funcs, 
              params@resource_dynamics,
              params@other_dynamics,
              params@other_encounter,
              params@other_mort,
              unique(params@gear_params$sel_func),
              kernel_fns)
    packages <- c("mizer", names(params@extensions))
    missing <- !sapply(packages, requireNamespace, quietly = TRUE)
    if (any(missing)) {
        stop("Some required extension packages are not installed: ",
             paste(missing, collapse = ", "))
    }
    custom <- sapply(funs, is_custom, packages = packages)
    if (any(custom)) {
        warning("Your model is using the functions ",
                paste(funs[custom], collapse = ", "),
                ". To share your model you need to share not only the ",
                "params object but also an R script or R Markdown file ",
                "defining these functions.")
    }
    saveRDS(params, file = file)
}

#' @rdname saveParams
#' @export
readParams <- function(file) {
    params <- validParams(readRDS(file))
    
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

is_custom <- function(name, packages) {
    !any(sapply(packages, 
                function(x, fun) exists(fun, where = paste0("package:", x),
                                        mode = "function"),
                fun = name))
}
