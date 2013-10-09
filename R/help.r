#' mizer: Multi-species size-based modelling in R
#'
#' The mizer package implements multi-species size-based modelling in R. It has been designed for modelling marine ecosystems.
#' 
#' Using \pkg{mizer} is relatively simple.  There are three main stages:
#'
#' Setting the model parameters. This is done by creating an object of class \code{MizerParams}. This includes model parameters such as the life history parameters of each species, and the range of the size spectrum.
#' Running a simulation. This is done by calling the \code{project()} method on the model parameters. This produces an object of \code{MizerSim} which contains the results of the simulation.
#' Exploring results. After a simulation has been run, the results can be explored using a range of plots and summaries.
#'
#' See the accompanying vignette for full details of the principles behind mizer and how the package can be used to peform size-based modelling.
#'
#' @import plyr ggplot2 grid methods
#' @importFrom reshape2 melt
#'
#' @docType package
#' @name mizer
#' @aliases mizer, mizer-package
NULL
