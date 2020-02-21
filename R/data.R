# NS_species_params ----
#' Example species parameter set based on the North Sea
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
#' \item{R_max}{Maximum recruitment}
#' \item{k_vb}{The von Bertalanffy k parameter}
#' }
#' @source{Blanchard et al.}
#' @docType data
#' @name NS_species_params
#' @examples 
#' \dontrun{
#' params <- MizerParams(NS_species_params)
#' sim = project(params)
#' plot(sim)
#' }
NULL

# NS_species_params_gears ----
#' Example species parameter set based on the North Sea with different gears
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
#' \item{R_max}{Maximum recruitment}
#' \item{k_vb}{The von Bertalanffy k parameter}
#' \item{gear}{Name of the fishing gear}
#' }
#' @source{Blanchard et al.}
#' @docType data
#' @name NS_species_params_gears
#' @examples 
#' \dontrun{
#' params <- MizerParams(NS_species_params_gears)
#' sim = project(params, effort = c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5))
#' plot(sim)
#' }
NULL

# inter ----
#' Example interaction matrix for the North Sea example
#'
#' The interaction coefficient between predators and preys in the North Sea.
#' @format A 12 x 12 matrix.
#' @source{Blanchard et al.}
#' @docType data
#' @name inter
NULL

# NS_params ----
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
#' @family example parameter objects
#' @examples 
#' \dontrun{
#' sim = project(NS_params, effort = c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5))
#' plot(sim)
#' }
NULL

# Benguela current ----
#' Example MizerParams object for the Benguela current
#' 
#' Created with purely size-based predation, i.e., no species-specific
#' interactions. Set up with three fishing gears targeting small, medium and
#' large species. Vulnerabilities are represented by changing the clearance rate
#' constant (gamma) between species. Calibrated to efforts
#' \code{effort = c(small = 0.13, medium = 0.05, large = 0.45)}.
#' @examples 
#' \dontrun{
#' sim = project(Benguela_params, effort = c(small = 0.13, medium = 0.05, large = 0.45))
#' plot(sim)
#' }
#' @format A MizerParams object
#' @source{N.S. Jacobsen, M. Burgess and K.H. Andersen (2017): 
#' Efficiency of fisheries is increasing at the ecosystem level. 
#' Fish and Fisheries 18(2) 199- 211. doi:10.1111/faf.12171.}
#' @docType data
#' @name Benguela_params
#' @family example parameter objects
NULL

# Baltic Sea ----
#' Example MizerParams object for the Central Baltic Sea
#' 
#' Created with purely size-based predation, i.e., no species-specific
#' interactions. Set up with three fishing gears targeting small, medium and
#' large species. Vulnerabilities are represented by changing the clearance rate
#' constant (gamma) between species. Calibrated to efforts
#' \code{effort = c(small = 0.3, medium = 0.3, large = 0.7)}.
#' @examples 
#' \dontrun{
#' sim = project(Baltic_params, effort = c(small = 0.3, medium = 0.3, large = 0.7))
#' plot(sim)
#' }
#' @format A MizerParams object
#' @source{N.S. Jacobsen, M. Burgess and K.H. Andersen (2017): 
#' Efficiency of fisheries is increasing at the ecosystem level. 
#' Fish and Fisheries 18(2) 199- 211. doi:10.1111/faf.12171.}
#' @docType data
#' @name Baltic_params
#' @family example parameter objects
NULL

# Barents Sea ----
#' Example MizerParams object for the Barents Sea
#' 
#' Created with purely size-based predation, i.e., no species-specific
#' interactions. Set up with three fishing gears targeting small, medium and
#' large species. Vulnerabilities are represented by changing the clearance rate
#' constant (gamma) between species. Calibrated to efforts
#' \code{effort = c(small = 1.1, medium = 0.5, large = 0.75)}.
#' @examples 
#' \dontrun{
#' sim = project(Barents_params, effort = c(small = 1.1, medium = 0.5, large = 0.75))
#' plot(sim)
#' }
#' @format A MizerParams object
#' @source{N.S. Jacobsen, M. Burgess and K.H. Andersen (2017): 
#' Efficiency of fisheries is increasing at the ecosystem level. 
#' Fish and Fisheries 18(2) 199- 211. doi:10.1111/faf.12171.}
#' @docType data
#' @name Barents_params
#' @family example parameter objects
NULL

# North-East US ----
#' Example MizerParams object for the North East US Continental Shelf (NEUSCS) with 24 species.
#' 
#' Created with purely size-based predation, i.e., no species-specific
#' interactions. Set up with three fishing gears targeting small, medium and
#' large species. Vulnerabilities are represented by changing the clearance rate
#' constant (gamma) between species. Calibrated to efforts
#' \code{effort = c(small = 0.4, medium = 0.3, large = 0.25)}.
#' @examples 
#' \dontrun{
#' sim = project(NEUSCS_params, effort = c(small = 0.4, medium = 0.3, large = 0.25))
#' plot(sim)
#' }
#' @format A MizerParams object
#' @source{N.S. Jacobsen, M. Burgess and K.H. Andersen (2017): 
#' Efficiency of fisheries is increasing at the ecosystem level. 
#' Fish and Fisheries 18(2) 199- 211. doi:10.1111/faf.12171.}
#' @docType data
#' @name NEUSCS_params
#' @family example parameter objects
NULL

# North Sea 10 species ----
#' Example MizerParams object for the North Sea with 10 species.
#' 
#' Created with purely size-based predation, i.e., no species-specific
#' interactions. Set up with three fishing gears targeting small, medium and
#' large species. Vulnerabilities are represented by changing the clearance rate
#' constant (gamma) between species. Calibrated to efforts
#' \code{effort = c(small = 0.6, medium = 0.6, large = 1.25)}.
#' @examples 
#' \dontrun{
#' sim = project(NorthSea_params, effort = c(small = 0.6, medium = 0.6, large = 1.25))
#' plot(sim)
#' }
#' @format A MizerParams object
#' @source{N.S. Jacobsen, M. Burgess and K.H. Andersen (2017): 
#' Efficiency of fisheries is increasing at the ecosystem level. 
#' Fish and Fisheries 18(2) 199- 211. doi:10.1111/faf.12171.}
#' @docType data
#' @name NorthSea_params
#' @family example parameter objects
NULL

