#' Save a MizerParams object
#' 
#' `This`saveParams()`` saves a MizerParams object including all functions that
#' are referenced in the MizerParams object. This can then be loaded back with
#' `loadParams()`.
#' 
#' 
#' 
#' @param params A MizerParams object
#' @param file The name of the file or a connection where the MizerParams object
#'   is saved to or read from.
#' @export
saveParams <- function(params, file) {
    funs <- c(params@rates_funcs, 
              params@resource_dynamics,
              params@other_dynamics,
              params@other_encounter,
              params@other_mort)
    packages <- c("mizer", names(params@extensions))
    missing <- !sapply(packages, requireNamespace, quietly = TRUE)
    if (any(missing)) {
        stop("Some required extension packages are not installed: ",
             paste(missing, collapse = ", "))
    }
    custom <- sapply(funs, is_custom, packages = packages)
    to_save <- as.vector(c("params", funs[custom]), mode = "character")
    save(list = to_save, file = file)
}

#' @rdname saveParams
#' @export
loadParams <- function(file) {
    load(file)
}

is_custom <- function(name, packages) {
    !any(sapply(packages, 
                function(x, fun) exists(fun, where = paste0("package:", x),
                                        mode = "function"),
                fun = name))
}