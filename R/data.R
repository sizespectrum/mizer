#' Example parameter set based on the North Sea
#' 
#' This data set is based on species in the North Sea (Blanchard et al.).
#' It is a data.frame that contains all the necessary information to be used by the
#' \code{\link{MizerParams}} constructor.
#' As there is no gear column, each species is assumed to be fished by a
#' separate gear.
#' 
#' @format A data frame with 12 rows and 7 columns. Each row is a species.
#' \describe{
#' \item{species}{Name of the species}
#' \item{w_inf}{The von Bertalanffy W_infinity parameter}
#' \item{w_mat}{Size at maturity}
#' \item{beta}{Size preference ratio}
#' \item{sigma}{Width of the size-preference}
#' \item{r_max}{Maximum recruitment}
#' \item{k_vb}{The von Bertalanffy k parameter}
#' }
#' @source{Blanchard et al.}
#' @docType data
#' @name NS_species_params
NULL
# "NS_species_params"

#' Example parameter set based on the North Sea with different gears
#'
#' This data set is based on species in the North Sea (Blanchard et al.).
#' It is similar to the data set \code{NS_species_params} except that
#' this one has an additional column specifying the fishing gear that
#' operates on each species. 
#' @format A data frame with 12 rows and 8 columns. Each row is a species.
#' \describe{
#' \item{species}{Name of the species}
#' \item{w_inf}{The von Bertalanffy W_infinity parameter}
#' \item{w_mat}{Size at maturity}
#' \item{beta}{Size preference ratio}
#' \item{sigma}{Width of the size-preference}
#' \item{r_max}{Maximum recruitment}
#' \item{k_vb}{The von Bertalanffy k parameter}
#' \item{gear}{Name of the fishing gear}
#' }
#' @source{Blanchard et al.}
#' @docType data
#' @name NS_species_params_gears
NULL
# "NS_species_params_gears"

#' Example interaction matrix for the North Sea example
#'
#' The interaction coefficient between predators and preys in the North Sea.
#' @format A 12 x 12 matrix.
#' @source{Blanchard et al.}
#' @docType data
#' @name inter
NULL
# "inter"

#' Example MizerParams object for the North Sea example
#'
#' A MizerParams object created from the \code{NS_species_params_gears} species
#' parameters and the \code{inter} interaction matrix together with an initial
#' condition corresponding to the steady state obtained from fishing with an
#' effort 
#' \code{effort = c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)}.
#' @format A MizerParams object
#' @source{Blanchard et al.}
#' @docType data
#' @name NS_params
NULL
# "NS_params"
