#' mizer: Multi-species size-based modelling in R
#'
#' The mizer package implements multi-species size-based modelling in R. It has 
#' been designed for modelling marine ecosystems.
#'
#' Using \pkg{mizer} is relatively simple.  There are three main stages: 
#' \enumerate{
#'
#' \item Setting the model parameters. This is done by creating an object of
#' class \linkS4class{MizerParams}. This includes model parameters such as the
#' life history parameters of each species, and the range of the size spectrum.
#' There are several setup functions that help to create a MizerParams objects
#' for particular types of models:
#' \itemize{
#'   \item \code{\link{set_community_model}}
#'   \item \code{\link{set_trait_model}}
#'   \item \code{\link{set_scaling_model}}
#'   \item \code{\link{set_multispecies_model}}
#' }
#' \item Running a simulation. This is done by calling the
#' \code{\link{project}} function with the model parameters. This produces an
#' object of \linkS4class{MizerSim} that contains the results of the simulation.
#'
#' \item Exploring results. After a simulation has been run, the results can be
#' explored using a range of \code{\link{plotting_functions}} and
#' \code{\link{summary_functions}}.
#' }
#'
#' See the mizer website and vignettes for full details of the principles behind
#' mizer and how the package can be used to perform size-based modelling.
#'
#' @import ggplot2 methods assertthat shiny dplyr
#' @importFrom reshape2 melt
#' @importFrom plotly ggplotly plotlyOutput renderPlotly
#' @importFrom stats fft mvfft lm pnorm runif complete.cases
#' @docType package
#' @name mizer
#' @aliases mizer-package
NULL
