#' Example parameters data sets 
#' 
#' There are two example data sets, \code{NS_species_params} and \code{NS_species_params_gears} based on species in the North Sea (Blanchard et el.)
#' They are both data frames and are identical except that one has an additional column specifying the fishing gear that operates on each species.
#' The data frames contain all the necessary information to be used by the \code{\link{MizerParams}} constructor.
#' The data set without a fishing gear column will set up a model in which each species is fished by a separate gear.
#'  
#' @rdname NS_species_params
#' @name NS_species_params
#' @aliases NS_species_params
#' @aliases NS_species_params_gears
#' @seealso \code{MizerParams}
#' @docType data
#' @references The North Sea paper (Blanchard et al)
NULL

