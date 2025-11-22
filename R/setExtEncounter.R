#' Set external encounter rate
#' 
#' @section Setting external encounter rate:
#' The external encounter rate is the rate at which a predator encounters
#' food that is not explicitly modelled. It is a rate with units mass/year.
#' 
#' The `ext_encounter` argument allows you to specify an external encounter rate
#' that depends on species and body size. You can see an example of this in
#' the Examples section of the help page for [setExtEncounter()].
#' 
#' @param params MizerParams
#' @param ext_encounter Optional. An array (species x size) holding the external
#'   encounter rate.  If not supplied, the external encounter rate is left
#'   unchanged. Initially is is set to 0.
#' @param ... Unused
#' 
#' @return `setExtEncounter()`: A MizerParams object with updated external encounter
#'   rate.
#' @export
#' @family functions for setting parameters
#' @examples
#' params <- newMultispeciesParams(NS_species_params)
#'
#' #### Setting allometric encounter rate #######################
#' 
#' # Set coefficient for each species. Here we choose 0.1 for each species
#' encounter_pre <- rep(0.1, nrow(species_params(params)))
#' 
#' # Multiply by power of size with exponent, here chosen to be 3/4
#' # The outer() function makes it an array species x size
#' allo_encounter <- outer(encounter_pre, w(params)^(3/4))
#' 
#' # Change the external encounter rate in the params object
#' ext_encounter(params) <- allo_encounter
setExtEncounter <- function(params, ext_encounter = NULL, ...) {
    UseMethod("setExtEncounter")
}
#' @export
setExtEncounter.MizerParams <- function(params, ext_encounter = NULL, ...) {
    assert_that(is(params, "MizerParams"))
    
    if (is.null(ext_encounter)) {
        ext_encounter <- params@ext_encounter
    }
    
    assert_that(is.array(ext_encounter),
                identical(dim(ext_encounter), dim(params@ext_encounter)))
    params@ext_encounter[] <- ext_encounter
    
    # Keep old comment if new comment is NULL
    if (!is.null(comment(ext_encounter))) {
        comment(params@ext_encounter) <- comment(ext_encounter)
    }
    
    params@time_modified <- lubridate::now()
    return(params)
}

#' @rdname setExtEncounter
#' @return `getExtEncounter()` or equivalently `ext_encounter()`: An array
#' (species x size) with the external encounter rate.
#' @export
getExtEncounter <- function(params) {
    UseMethod("getExtEncounter")
}
#' @export
getExtEncounter.MizerParams <- function(params) {
    params@ext_encounter
}

#' @rdname setExtEncounter
#' @export
ext_encounter <- function(params) {
    UseMethod("ext_encounter")
}
#' @export
ext_encounter.MizerParams <- function(params) {
    params@ext_encounter
}

#' @rdname setExtEncounter
#' @param value ext_encounter
#' @export
`ext_encounter<-` <- function(params, value) {
    UseMethod("ext_encounter<-")
}
#' @export
`ext_encounter<-.MizerParams` <- function(params, value) {
    setExtEncounter(params, ext_encounter = value)
}