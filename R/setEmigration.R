#' Set emigration rate
#'
#' @section Setting emigration rate:
#' The emigration rate is the rate at which the abundance density is decreased over
#' time beyond the decrease described by the per-capita mortality.
#'
#' The `emigration` argument allows you to specify an emigration rate that depends on
#' species and body size.
#'
#' @param params MizerParams
#' @param emigration Optional. An array (species x size) holding the emigration rate. If
#'   not supplied, the emigration rate is left unchanged. Initially it is set to 0.
#' @param ... Unused
#'
#' @return `setEmigration()`: A MizerParams object with updated emigration rate.
#' @export
#' @family functions for setting parameters
setEmigration <- function(params, emigration = NULL, ...) {
    assert_that(is(params, "MizerParams"))

    if (is.null(emigration)) {
        emigration <- params@emigration
    }

    assert_that(is.array(emigration),
                identical(dim(emigration), dim(params@emigration)))
    params@emigration[] <- emigration

    # Keep old comment if new comment is NULL
    if (!is.null(comment(emigration))) {
        comment(params@emigration) <- comment(emigration)
    }

    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setEmigration
#' @return `emigration()`: An array (species x size) with the emigration rate.
#' @export
emigration <- function(params) {
    params@emigration
}

#' @rdname setEmigration
#' @param value emigration
#' @export
`emigration<-` <- function(params, value) {
    setEmigration(params, emigration = value)
} 