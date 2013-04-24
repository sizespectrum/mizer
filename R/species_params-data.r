#' Example parameters data sets 
#' 
#' There are two example data sets, \code{species_params_gears} and \code{species_params_nogears}.
#' They are both data frames and are identical except that one has an additional column specifying the fishing gear that operates on each species.
#' The data frames contain all the necessary information to be used by the \code{\link{MizerParams}} constructor.
#' The data set without a fishing gear column will set up a model in which each species is fished by a separate gear.
#'  
#' @rdname species_params
#' @name species_params
#' @aliases species_params_nogears
#' @aliases species_params_gears
#' @seealso \code{MizerParams}
#' @docType data
#' @references The North Sea paper (Blanchard et al)
#' @keywords species_params_nogears
NULL

