# Documentation for methods of base generics
#
# Copyright 2026 Gustav Delius.
# Distributed under the GPL 3 or later.

# print ----
#' Print mizer objects
#'
#' Mizer supplies `print()` methods for the array-like objects returned by many
#' rate and summary functions. These methods print a compact, readable overview
#' instead of the full matrix or array. The printed output reports the value
#' name, dimensions, units if known, and per-species minimum, mean and maximum
#' values.
#'
#' For full numeric access, use the object itself as an ordinary matrix or array
#' or convert it to a long data frame with [as.data.frame()].
#'
#' @param x The object to print.
#' @param ... Further arguments. They are currently ignored by the mizer
#'   methods.
#'
#' @return The printed object, invisibly.
#'
#' @aliases print.ArraySpeciesBySize print.ArrayTimeBySpecies print.ArrayTimeBySpeciesBySize print.summary.ArraySpeciesBySize print.summary.ArrayTimeBySpecies print.summary.ArrayTimeBySpeciesBySize
#' @seealso [summary()], [as.data.frame()], [plot()],
#'   [ArraySpeciesBySize()], [ArrayTimeBySpecies()],
#'   [ArrayTimeBySpeciesBySize()]
#'
#' @examples
#' \donttest{
#' enc <- getEncounter(NS_params)
#' print(enc)
#'
#' biomass <- getBiomass(NS_sim)
#' print(biomass)
#' }
#' @name print
NULL

# summary ----
#' Summarise mizer objects
#'
#' Mizer provides `summary()` methods for model objects and for the specialised
#' array classes returned by many mizer functions.
#'
#' For a [MizerParams()] object, `summary()` prints the model metadata, size
#' grids, selected species parameters and fishing gear details. For a
#' [MizerSim()] object, it first prints the parameter summary and then reports
#' the simulated time period and output interval.
#'
#' For [ArraySpeciesBySize()], [ArrayTimeBySpecies()] and
#' [ArrayTimeBySpeciesBySize()] objects, `summary()` returns a small list with
#' the value name, units, dimensions and a per-species data frame containing
#' minimum, mean and maximum values. Printing that summary object gives the same
#' compact table in a human-readable form.
#'
#' @param object The object to summarise.
#' @param ... Further arguments. They are currently ignored by the mizer
#'   methods.
#'
#' @return
#' For [MizerParams()] and [MizerSim()], the object is returned invisibly.
#' For array objects, a list of class `summary.ArraySpeciesBySize`,
#' `summary.ArrayTimeBySpecies` or `summary.ArrayTimeBySpeciesBySize`.
#'
#' @aliases summary.ArraySpeciesBySize summary.ArrayTimeBySpecies summary.ArrayTimeBySpeciesBySize summary.MizerSim summary.MizerParams
#' @seealso [print()], [as.data.frame()], [str()], [MizerParams()], [MizerSim()],
#'   [ArraySpeciesBySize()], [ArrayTimeBySpecies()],
#'   [ArrayTimeBySpeciesBySize()]
#'
#' @examples
#' \donttest{
#' summary(NS_params)
#' summary(NS_sim)
#' summary(getEncounter(NS_params))
#' summary(getFMort(NS_sim))
#' }
#' @name summary
NULL

# as.data.frame ----
#' Convert mizer arrays to data frames
#'
#' The `as.data.frame()` methods for mizer array classes turn matrix- and
#' array-like results into tidy long-form data frames, with one row per
#' observed combination of species, size and/or time. The numeric result is
#' always stored in a column called `value`.
#'
#' The returned columns are:
#' \itemize{
#'   \item `ArraySpeciesBySize`: `w`, `value`, `Species`.
#'   \item `ArrayTimeBySpecies`: `time`, `value`, `Species`.
#'   \item `ArrayTimeBySpeciesBySize`: `time`, `Species`, `w`, `value`.
#' }
#'
#' If the original object has non-numeric or missing dimension names, sequential
#' indices are used for the `time` or `w` columns. Species names are taken from
#' the row, column or dimension names of the original object.
#'
#' @param x An `ArraySpeciesBySize`, `ArrayTimeBySpecies` or
#'   `ArrayTimeBySpeciesBySize` object.
#' @param row.names Optional and included only for compatibility with the base
#'   generic. `NULL` or a character vector giving the row names for the
#'   data frame.
#' @param optional Optional and included only for compatibility with the base
#'   generic. A logical value. If `TRUE`, setting row names and converting
#'   column names (to syntactic names) is optional.
#' @param ... Further arguments. They are currently ignored by the mizer
#'   methods.
#'
#' @return A data frame in long format.
#'
#' @aliases as.data.frame.ArraySpeciesBySize as.data.frame.ArrayTimeBySpecies as.data.frame.ArrayTimeBySpeciesBySize
#' @seealso [print()], [summary()], [str()], [plot()], [ArraySpeciesBySize()],
#'   [ArrayTimeBySpecies()], [ArrayTimeBySpeciesBySize()]
#'
#' @examples
#' \donttest{
#' enc <- getEncounter(NS_params)
#' head(as.data.frame(enc))
#'
#' biomass <- getBiomass(NS_sim)
#' head(as.data.frame(biomass))
#' }
#' @name as.data.frame
NULL

# str ----
#' Display the structure of mizer objects
#'
#' Mizer provides `str()` methods for [MizerParams()] and [MizerSim()] objects,
#' as well as [ArraySpeciesBySize()], [ArrayTimeBySpecies()] and [ArrayTimeBySpeciesBySize()]
#' objects. These methods produce a clean, compact overview of the object's structure
#' without polluting the console with large amounts of internal data.
#'
#' @param object The object to display the structure of.
#' @param max.level Maximum level of nesting to print. Defaults to `NA` (no limit).
#' @param ... Further arguments. They are passed to the default `str()` method.
#'
#' @return `NULL`, invisibly.
#'
#' @aliases str.ArraySpeciesBySize str.ArrayTimeBySpecies str.ArrayTimeBySpeciesBySize str.MizerSim str.MizerParams
#' @seealso [print()], [as.data.frame()], [summary()], [plot()], [MizerParams()], [MizerSim()], [ArraySpeciesBySize()],
#'   [ArrayTimeBySpecies()], [ArrayTimeBySpeciesBySize()]
#'
#' @examples
#' \donttest{
#' str(NS_params)
#' str(NS_sim)
#' str(getEncounter(NS_params))
#' }
#' @name str
NULL

