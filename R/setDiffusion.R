#' Set diffusion rate
#'
#' @section Setting diffusion rate:
#' The diffusion rate governs the redistribution of abundance density over body size
#' beyond deterministic growth dynamics.
#'
#' The `diffusion` argument allows you to specify a diffusion rate that depends on
#' species and body size.
#'
#' @param params MizerParams
#' @param diffusion Optional. An array (species x size) holding the diffusion rate. If
#'   not supplied, the diffusion rate is left unchanged. Initially it is set to 0.
#' @param reset `r lifecycle::badge("experimental")`
#'   If set to TRUE, then the diffusion rate will be reset to the value
#'   calculated from the species parameters, even if it was previously
#'   overwritten with a custom value. If set to FALSE (default) then a
#'   recalculation from the species parameters will take place only if no
#'   custom value has been set.
#' @param ... Unused
#'
#' @return `setDiffusion()`: A MizerParams object with updated diffusion rate.
#' @export
#' @family functions for setting parameters
setDiffusion <- function(params, diffusion = NULL, reset = FALSE, ...) {
    UseMethod("setDiffusion")
}
#' @rdname setDiffusion
#' @export
setDiffusion.MizerParams <- function(params, diffusion = NULL, reset = FALSE, ...) {
    assert_that(is(params, "MizerParams"))

    if (is.null(diffusion)) {
        diffusion <- params@diffusion
    }

    assert_that(is.array(diffusion),
                identical(dim(diffusion), dim(params@diffusion)))
    params@diffusion[] <- diffusion

    # Keep old comment if new comment is NULL
    if (!is.null(comment(diffusion))) {
        comment(params@diffusion) <- comment(diffusion)
    }

    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setDiffusion
#' @return `diffusion()`: An array (species x size) with the diffusion rate.
#' @export
diffusion <- function(params) {
    UseMethod("diffusion")
}
#' @export
diffusion.MizerParams <- function(params) {
    params@diffusion
}

#' @rdname setDiffusion
#' @param value diffusion
#' @export
`diffusion<-` <- function(params, value) {
    UseMethod("diffusion<-")
}
#' @export
`diffusion<-.MizerParams` <- function(params, value) {
    setDiffusion(params, diffusion = value)
}
