# NS_species_params ----
#' Example species parameter set based on the North Sea
#' 
#' This data set is based on species in the North Sea (Blanchard et al.). It is
#' a data.frame that contains all the necessary information to be used by the
#' [MizerParams()] constructor. As there is no gear column, each species is
#' assumed to be fished by a separate gear.
#' 
#' @format A data frame with 12 rows and 7 columns. Each row is a species.
#' \describe{
#' \item{species}{Name of the species}
#' \item{w_max}{The size at which the population invests 100% of its income into
#'   reproduction so that all growth stops.}
#' \item{w_mat}{Size at maturity}
#' \item{beta}{Size preference ratio}
#' \item{sigma}{Width of the size-preference}
#' \item{R_max}{Maximum reproduction rate}
#' \item{k_vb}{The von Bertalanffy k parameter}
#' \item{w_inf}{The von Bertalanffy asymptotic size}
#' }
#' @source{Blanchard et al.}
#' @docType data
#' @name NS_species_params
#' @examples
#' params <- MizerParams(NS_species_params)
"NS_species_params"

# NS_species_params_gears ----
#' Example species parameter set based on the North Sea with different gears
#'
#' This data set is based on species in the North Sea (Blanchard et al.).
#' It is similar to the data set `NS_species_params` except that
#' this one has an additional column specifying the fishing gear that
#' operates on each species. 
#' @format A data frame with 12 rows and 8 columns. Each row is a species.
#' \describe{
#' \item{species}{Name of the species}
#' \item{w_max}{The size at which the population invests 100% of its income into
#'   reproduction so that all growth stops.}
#' \item{w_mat}{Size at maturity}
#' \item{beta}{Size preference ratio}
#' \item{sigma}{Width of the size-preference}
#' \item{R_max}{Maximum reproduction rate}
#' \item{k_vb}{The von Bertalanffy k parameter}
#' \item{w_inf}{The von Bertalanffy asymptotic size}
#' \item{gear}{Name of the fishing gear}
#' }
#' @source{Blanchard et al.}
#' @docType data
#' @name NS_species_params_gears
#' @examples
#' params <- MizerParams(NS_species_params_gears)
"NS_species_params_gears"

# NS_interaction ----
#' Example interaction matrix for the North Sea example
#'
#' The interaction coefficient between predator and prey species
#' in the North Sea.
#' @format A 12 x 12 matrix.
#' @source{Blanchard et al.}
#' @docType data
#' @name NS_interaction
#' @examples
#' params <- MizerParams(NS_species_params_gears,
#'                       interaction = NS_interaction)
"NS_interaction"

#' Alias for `NS_interaction`
#' 
#' `r lifecycle::badge("deprecated")`
#' An alias provided for backward compatibility with mizer version <= 2.3
#' @format A 12 x 12 matrix.
#' @source{Blanchard et al.}
#' @docType data
#' @concept deprecated
"inter"

# NS_params ----
#' Example MizerParams object for the North Sea example
#'
#' A MizerParams object created from the `NS_species_params_gears` species
#' parameters and the `inter` interaction matrix together with an initial
#' condition corresponding to the steady state obtained from fishing with an
#' effort 
#' \code{effort = c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)}.
#' @format A MizerParams object
#' @source{Blanchard et al.}
#' @docType data
#' @name NS_params
#' @family example parameter objects
#' @examples 
#' \donttest{
#' sim = project(NS_params, effort = c(Industrial = 0, Pelagic = 1, 
#'                                     Beam = 0.5, Otter = 0.5))
#' plot(sim)
#' }
"NS_params"

# NS_sim ----
#' Example MizerSim object for the North Sea example
#'
#' A MizerSim object containing a simulation with historical fishing
#' mortalities from the North Sea, as created in the tutorial
#' "A Multi-Species Model of the North Sea".
#' @format A MizerSim object
#' @source \url{https://sizespectrum.org/mizer/articles/a_multispecies_model_of_the_north_sea.html}
#' @docType data
#' @name NS_sim
#' @family example parameter objects
#' @examples
#' plotBiomass(NS_sim)
"NS_sim"
