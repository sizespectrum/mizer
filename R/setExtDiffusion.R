#' Set external diffusion rate
#'
#' @section Setting external diffusion rate:
#' The external diffusion rate allows you to impose additional diffusion of
#' abundance density over body size beyond the predation-driven diffusion
#' computed by [mizerDiffusion()].
#'
#' The `ext_diffusion` argument allows you to specify a diffusion rate that
#' depends on species and body size.
#'
#' @param params MizerParams
#' @param ext_diffusion Optional. An array (species x size) holding the
#'   external diffusion rate. If not supplied, the external diffusion rate is
#'   left unchanged. Initially it is set to 0.
#' @param reset Unused. Included for interface consistency with other setter
#'   functions.
#' @param ... Unused
#'
#' @return `setExtDiffusion()`: A MizerParams object with updated external
#'   diffusion rate.
#' @export
#' @family functions for setting parameters
setExtDiffusion <- function(params, ext_diffusion = NULL, reset = FALSE, ...) {
    UseMethod("setExtDiffusion")
}
#' @rdname setExtDiffusion
#' @export
setExtDiffusion.MizerParams <- function(params, ext_diffusion = NULL,
                                        reset = FALSE, ...) {
    if (is.null(ext_diffusion)) {
        ext_diffusion <- params@ext_diffusion
    }

    assert_that(is.array(ext_diffusion),
                identical(dim(ext_diffusion), dim(params@ext_diffusion)))
    params@ext_diffusion[] <- ext_diffusion

    # Keep old comment if new comment is NULL
    if (!is.null(comment(ext_diffusion))) {
        comment(params@ext_diffusion) <- comment(ext_diffusion)
    }

    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setExtDiffusion
#' @return `ext_diffusion()`: An `ArraySpeciesBySize` object (species x size)
#'   with the external diffusion rate.
#' @export
ext_diffusion <- function(params) {
    UseMethod("ext_diffusion")
}
#' @export
ext_diffusion.MizerParams <- function(params) {
    ArraySpeciesBySize(params@ext_diffusion,
                       value_name = "Diffusion rate",
                       params = params)
}

#' @rdname setExtDiffusion
#' @param value ext_diffusion
#' @export
`ext_diffusion<-` <- function(params, value) {
    UseMethod("ext_diffusion<-")
}
#' @export
`ext_diffusion<-.MizerParams` <- function(params, value) {
    setExtDiffusion(params, ext_diffusion = value)
}
