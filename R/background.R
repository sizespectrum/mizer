#' Designate species as background species
#'
#' Marks the specified set of species as background species by setting the
#' `is_background` column in their species parameters to `TRUE`. Background
#' species are handled differently in plots (displayed in grey) and their
#' abundances can be automatically adjusted to keep the community close to the
#' Sheldon spectrum (see `adjustBackgroundSpecies()` in the mizerExperimental
#' package).
#'
#' @param object An object of class \linkS4class{MizerParams} or
#'   \linkS4class{MizerSim}.
#' @inheritParams valid_species_arg
#'
#' @return An object of the same class as the `object` argument
#' @export
#' @seealso [removeBackgroundSpecies()]
#' @examples
#' params <- markBackground(NS_params,
#'                          species = c("Sprat", "Sandeel", "N.pout"))
#' any(species_params(params)$is_background)
markBackground <- function(object, species = NULL) {
    if (is(object, "MizerSim")) {
        species <- valid_species_arg(object, species)
        idx <- dimnames(object@params@initial_n)$sp %in% species
        object@params@species_params$is_background[idx] <- TRUE
    } else if (is(object, "MizerParams")) {
        species <- valid_species_arg(object, species)
        idx <- dimnames(object@initial_n)$sp %in% species
        object@species_params$is_background[idx] <- TRUE
    } else {
        stop("The `object` argument must be of type MizerParams or MizerSim.")
    }
    return(object)
}


#' Remove all background species
#'
#' Removes all species that have been marked as background species with
#' [markBackground()].
#' 
#' This is just a shorthand for
#' `removeSpecies(params, species_params(params)$is_background)`
#'
#' @param params A \linkS4class{MizerParams} object
#' @return A \linkS4class{MizerParams} object with background species removed
#' @export
#' @seealso [markBackground()]
#' @examples
#' params <- markBackground(NS_params,
#'                          species = c("Sprat", "Sandeel", "N.pout"))
#' params <- removeBackgroundSpecies(params)
#' species_params(params)$species
removeBackgroundSpecies <- function(params) {
    removeSpecies(params, params@species_params$is_background)
}
