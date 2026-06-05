#' @param object A \linkS4class{MizerParams} or \linkS4class{MizerSim} object.
#' @param ... Additional arguments that depend on the class of `object`.
#'
#'   **For a \linkS4class{MizerParams} object:**
#'   \describe{
#'     \item{`n`}{A matrix of species abundances (species x size). Defaults to
#'       the initial abundances stored in `object`.}
#'     \item{`n_pp`}{A vector of the resource abundance by size. Defaults to the
#'       initial resource abundance stored in `object`.}
#'     \item{`n_other`}{A named list of the abundances of other dynamical
#'       components. Defaults to the initial values stored in `object`.}
#'     \item{`t`}{The time for which to do the calculation. Defaults to 0.}
#'   }
#'
#'   **For a \linkS4class{MizerSim} object:**
#'   \describe{
#'     \item{`time_range`}{The time range over which to return the rates. Either
#'       a vector of values, a vector of min and max time, or a single value.
#'       Defaults to the whole time range of the simulation.}
#'   }
