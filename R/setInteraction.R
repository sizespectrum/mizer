#' Set species interaction matrix
#'
#' @section Setting interactions:
#' 
#' The interaction matrix \eqn{\theta_{ij}} describes the interaction of each
#' pair of species in the model. This can be viewed as a proxy for spatial
#' interaction e.g. to model predator-prey interaction that is not size based.
#' The values in the interaction matrix are used to scale the encountered food
#' and predation mortality (see on the website
#' [the section on predator-prey encounter rate](https://sizespectrum.org/mizer/articles/model_description.html#sec:pref)
#' and on [predation mortality](https://sizespectrum.org/mizer/articles/model_description.html#mortality)).
#'
#' It is used when calculating the food encounter rate in [getEncounter()] and
#' the predation mortality rate in [getPredMort()]. Its entries are
#' dimensionless numbers.The values are between 0 (species do not overlap and
#' therefore do not interact with each other) to 1 (species overlap perfectly).
#' If all the values in the interaction matrix are set to 1 then predator-prey
#' interactions are determined entirely by size-preference.
#' 
#' This function checks that the supplied interaction matrix is valid and then
#' stores it in the `interaction` slot of the params object before returning
#' that object.
#'
#' The order of the columns and rows of the `interaction` argument should be the
#' same as the order in the species params data frame in the `params` object. If
#' you supply a named array then the function will check the order and warn if
#' it is different. One way of creating your own interaction matrix is to enter
#' the data using a spreadsheet program and saving it as a .csv file. The data
#' can be read into R using the command `read.csv()`.
#'
#' The interaction of the species with the resource are set via a column
#' `interaction_p` in the `species_params` data frame. Again the entries have to
#' be numbers between 0 and 1. By default this column is set to all 1s.
#'
#' @param params MizerParams object
#' @param interaction Optional interaction matrix of the species (predator
#'   species x prey species). Entries should be numbers between 0 and 1. By
#'   default all entries are 1. See "Setting interactions" section below.
#'
#' @return MizerParams object with updated interaction matrix. Because of the
#'   way the R language works, `setInteraction()` does not make the changes to
#'   the params object that you pass to it but instead returns a new params
#'   object. So to affect the change you call the function in the form
#'   `params <- setInteraction(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- newTraitParams()
#' interaction <- params@interaction
#' interaction[1, 3] <- 0
#' params <- setInteraction(params, interaction)
#' }
setInteraction <- function(params,
                           interaction = NULL) {
    assert_that(is(params, "MizerParams"))
    if (is.null(interaction)) {
        interaction <- params@interaction
    }
    if (!is.matrix(interaction)) {
        interaction <- as.matrix(interaction)
    }
    # Check dims of interaction argument
    if (!identical(dim(params@interaction), dim(interaction))) {
        stop("interaction matrix is not of the right dimensions. Must be ",
             "number of species x number of species.")
    }
    # Check that all values of interaction matrix are 0 - 1.
    if (!all((interaction >= 0) & (interaction <= 1))) {
        stop("Values in the interaction matrix must be between 0 and 1")
    }
    # In case user has supplied names to interaction matrix, check them.
    if (!is.null(dimnames(interaction))) {
        if (!is.null(names(dimnames(interaction)))) {
            if (!identical(names(dimnames(interaction)),
                           names(dimnames(params@interaction)))) {
                message("Note: Your interaction matrix has dimensions called: `",
                        toString(names(dimnames(interaction))),
                        "`. I expected 'predator, prey'. ", 
                        "I will now ignore your names.")
            }
        }
        names(dimnames(interaction)) <- names(dimnames(params@interaction))
        if (!identical(dimnames(params@interaction),
                       dimnames(interaction))) {
            message("Note: Dimnames of interaction matrix do not match the ",
                    "order of species names in the species data.frame. I am ",
                    "now ignoring your dimnames so your interaction matrix ",
                    "may be in the wrong order.")
        }
    }
    params@interaction[] <- interaction
    
    # Check the interaction_p column in species_params
    message <- "Note: No interaction_p column in species data frame so assuming all species feed on resource."
    species_params <- set_species_param_default(params@species_params,
                                                "interaction_p", 1,
                                                message = message)
    # Check that all values of interaction vector are 0 - 1.
    if (!all((species_params$interaction_p >= 0) & 
             (species_params$interaction_p <= 1))) {
        stop("Values in the resource interaction vector should be between 0 and 1")
    }
    params@species_params$interaction_p <- species_params$interaction_p
    
    return(params)
}

#' @rdname setInteraction
#' @export
getInteraction <- function(params) {
    params@interaction
}
