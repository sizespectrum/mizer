#' Set Beverton-Holt reproduction without changing the steady state
#'
#' Takes a MizerParams object `params` with arbitrary density dependence in
#' reproduction and
#' returns a MizerParams object with Beverton-Holt density dependence while
#' preserving the current effective recruitment at the initial state. Hence if
#' you have tuned your `params` object to describe a
#' particular steady state, then setting the Beverton-Holt density dependence
#' with this function will leave you with the exact same steady state. By
#' specifying one of the parameters `erepro`, `R_max` or `reproduction_level`
#' you pick the desired reproduction curve. More details of these parameters are
#' provided below.
#'
#' With Beverton-Holt density dependence the relation between the energy
#' invested into reproduction and the number of eggs hatched is determined
#' by two parameters: the reproductive efficiency `erepro` and the maximum
#' reproduction rate `R_max`.
#'
#' If no maximum is imposed on the reproduction rate
#' (\eqn{R_{max} = \infty}{R_max = Inf}) then the resulting density-independent
#' reproduction rate \eqn{R_{di}}{R_di} is proportional
#' to the total rate \eqn{E_R} at which energy is invested into reproduction,
#' \deqn{R_{di} = \frac{\rm{erepro}}{2 w_{min}} E_R,}{
#'   R_di = (erepro/(2 w_min)) E_R,}
#' where the proportionality factor is given by the reproductive efficiency
#' `erepro` divided by the egg size `w_min` to convert energy to egg number and
#' divided by 2 to account for the two sexes.
#'
#' Imposing a finite maximum reproduction rate \eqn{R_{max}}{R_max} leads to a
#' non-linear relationship between energy invested and eggs hatched. This
#' density-dependent reproduction rate \eqn{R_{dd}}{R_dd} is given as
#' \deqn{R_{dd} = R_{di}
#' \frac{R_{max}}{R_{di} + R_{max}}.}{R_dd = R_di R_max/(R_di + R_max).}
#'
#' (All quantities in the above equations are species-specific but we dropped
#' the species index for simplicity.)
#'
#' The following plot illustrates the Beverton-Holt density dependence in the
#' reproduction rate for two different choices of parameters.
#'
#' ```{r Beverton-Holt-plot, echo = FALSE, fig.height = 4.5, fig.width = 9}
#' erepro <- 4
#' R_max <- 1
#' E_R <- seq(0, 2, by = 0.05)
#' R_di = erepro * E_R
#' R_dd <- R_di * R_max / (R_di + R_max)
#' df <- melt(data.frame(E_R, R_dd, R_di, R_max), id.vars = "E_R")
#' df <- df[df$value < 1.6, ]
#' df$dd <- "Low"
#'
#' erepro <- 1.5
#' R_max <- 3/2
#' R_di = erepro * E_R
#' R_dd <- R_di * R_max / (R_di + R_max)
#' df2 <- melt(data.frame(E_R, R_dd, R_di, R_max), id.vars = "E_R")
#' df2 <- df2[df2$value < 1.6, ]
#' df2$dd <- "High"
#'
#' ggplot(rbind(df, df2)) +
#'     geom_line(aes(x = E_R, y = value, linetype = variable,
#'                   colour = dd, linewidth = dd)) +
#'     annotate("point", x = 5/4, y = 5/6, size = 2) +
#'     labs(linetype = "", linewidth = "R_max", colour = "R_max") +
#'     scale_y_continuous(name = "Reproduction rate [eggs/year]",
#'                        breaks = c(5/6), labels = c("R_dd")) +
#'     scale_x_continuous(name = "Energy invested [g/year]",
#'                        breaks = c(5/4), labels = c("E_R")) +
#'     scale_linewidth_manual(values = c(1, 0.5)) +
#'     scale_colour_manual(values = c("black", "blue")) +
#'     scale_linetype_manual(values = c("solid", "dashed", "dotted"))
#' ```
#'
#' This plot shows that a given energy \eqn{E_R} invested into reproduction can
#' lead to the same reproduction rate \eqn{R_{dd}}{R_dd} with different choices
#' of the parameters `R_max` and `erepro`. `R_max` determines the asymptote of
#' the curve and `erepro` its initial slope. A higher `R_max` coupled with a
#' lower `erepro` (black curves) can give the same value as a lower `R_max`
#' coupled with a higher `erepro` (blue curves).
#'
#' For the given initial state in the MizerParams object `params` one can
#' calculate the energy \eqn{E_R} invested into reproduction and the current
#' effective recruitment \eqn{R_{dd}}{R_dd}. These two values determine the
#' location of the black dot in the above graph. You then only need one
#' parameter to select one curve from the family of Beverton-Holt curves going
#' through that point. This parameter can be `erepro` or `R_max`. Instead of
#' `R_max` you can alternatively specify the `reproduction_level` which is the
#' ratio between the density-dependent reproduction rate \eqn{R_{dd}}{R_dd} and
#' the maximal reproduction rate  \eqn{R_{max}}{R_max}.
#'
#' If you do not provide a value for any of the reproduction parameter
#' arguments, then `erepro` will be set to the value it has in the current
#' species parameter data frame. If you do provide one of the reproduction
#' parameters, this can be either a vector with one value for each
#' species, or a named vector where the names determine which species are
#' affected, or a single unnamed value that is then used for all species. Any
#' species for which the given value is `NA` will remain unaffected.
#'
#' The values for `R_max` must be larger than \eqn{R_{dd}}{R_dd} and can range
#' up to `Inf`. An infeasible requested value is an error.
#'
#' The values for the `reproduction_level` must be non-negative and
#' less than 1. The values for `erepro` must be large enough to allow the
#' preserved recruitment rate. A smaller requested value is an error. The values
#' for `erepro` should also be smaller than 1 to be physiologically sensible,
#' but this is not enforced by the function.
#'
#' A species that is absent from the initial state, with zero abundance at all
#' sizes, has no reproduction whatever its parameters, so every Beverton-Holt
#' curve keeps it at its steady state. For such a species a given `erepro` or
#' `R_max` is used as it is and the other parameter is left unchanged. Its
#' reproduction level is zero whatever its `R_max`, so a non-zero
#' `reproduction_level` cannot be set for it and is rejected.
#'
#' As can be seen in the graph above, choosing a lower value for `R_max` or a
#' higher value for `erepro` means that near the steady state the reproduction
#' will be less sensitive to a change in the energy invested into reproduction
#' and hence less sensitive to changes in the spawning stock biomass or its
#' energy income. As a result the species will also be less sensitive to
#' fishing, leading to a higher F_MSY.
#'
#' @param params A MizerParams object
#' @param erepro Reproductive efficiency for each species. See details.
#' @param R_max Maximum reproduction rate. See details.
#' @param reproduction_level Sets `R_max` so that the reproduction rate at
#'   the initial state is `R_max * reproduction_level`.
#' @param info_level Controls the amount of information messages and warnings
#'   that are shown. Higher levels lead to more messages, `info_level = 0`
#'   gives silence. The default is taken from the `mizer_info_level` option,
#'   see [default_info_level()].
#' @param ... Unused
#'   \itemize{
#'     \item `R_factor`: Legacy alternative for specifying
#'       `reproduction_level = 1 / R_factor`.
#'   }
#'
#' @return A MizerParams object
#' @examples
#' params <- NS_params
#' species_params(params)$erepro
#' # Increase efficiency while preserving recruitment
#' params <- setBevertonHolt(params,
#'     erepro = 2 * species_params(params)$erepro
#' )
#' t(species_params(params)[, c("erepro", "R_max")])
#' # Setting erepro for some species
#' params <- setBevertonHolt(params, reproduction_level = c(Gurnard = 0.6))
#' t(species_params(params)[, c("erepro", "R_max")])
#' # Setting R_max
#' R_max <- 2 * getRDD(params)
#' params <- setBevertonHolt(NS_params, R_max = R_max)
#' t(species_params(params)[, c("erepro", "R_max")])
#' # Setting reproduction_level
#' params <- setBevertonHolt(params, reproduction_level = 0.3)
#' t(species_params(params)[, c("erepro", "R_max")])
#' @export
setBevertonHolt <- function(params, erepro,
                            # nolint next: object_name_linter.
                            R_max, reproduction_level,
                            info_level = default_info_level(), ...) {
    UseMethod("setBevertonHolt")
}
#' @export
setBevertonHolt.MizerParams <- function(params, erepro,
                                        # nolint next: object_name_linter.
                                        R_max, reproduction_level,
                                        info_level = default_info_level(),
                                        ...) {
    with_info_level(info_level = info_level, {
        params <- validParams(params)
        controls <- list(...)
        unsupported <- setdiff(names(controls), "R_factor")
        if (length(unsupported)) {
            stop("Unsupported Beverton-Holt arguments: ",
                paste(unsupported, collapse = ", "), ".",
                call. = FALSE
            )
        }
        if (!missing(erepro)) controls$erepro <- erepro
        if (!missing(R_max)) controls$R_max <- R_max
        if (!missing(reproduction_level)) {
            controls$reproduction_level <- reproduction_level
        }
        if (length(controls) > 1L) {
            stop("You should only provide `params` and one other argument.",
                call. = FALSE
            )
        }
        if (!length(controls)) controls$erepro <- params@species_params$erepro
        field <- names(controls)[[1L]]
        value <- reproduction_species_values(params, controls[[1L]], field,
            allow_infinite = field %in% c(
                "R_max",
                "R_factor"
            )
        )
        if (field == "R_factor") {
            if (any(value <= 1, na.rm = TRUE)) {
                stop("The R_factor must be greater than 1.", call. = FALSE)
            }
            value <- 1 / value
            field <- "reproduction_level"
        }
        selected <- names(value)[!is.na(value)]
        if (!length(selected)) {
            return(params)
        }
        if (field == "reproduction_level" &&
                any(value < 0 | value >= 1, na.rm = TRUE)) {
            stop("The reproduction level must be smaller than 1 ",
                "and non-negative.",
                call. = FALSE
            )
        }
        switching <- reproduction_law_name(params) != "BevertonHoltRDD"
        if (switching && length(selected) != nrow(params@species_params)) {
            stop("Switching the reproduction law is model-wide; supply ",
                "Beverton-Holt controls for every species.",
                call. = FALSE
            )
        }
        rates <- tuning_rate_context(params, params@initial_effort, 0)
        targets <- stats::setNames(
            as.numeric(rates$rdd),
            params@species_params$species
        )[selected]
        if (switching) {
            params <- set_species_param_default(params, "R_max", Inf)
            params <- setRateFunction(params, "RDD", "BevertonHoltRDD")
        }
        args <- c(
            list(params, rdd_target = targets),
            stats::setNames(list(value), field),
            list(info_level = info_level)
        )
        do.call(balanceRDD, args)
    })
}
#' @rdname setBevertonHolt
#' @return `reproduction_level()`: A named vector with the reproduction level
#'   for each species.
#' @export
#' @examples
#'
#' # Inspecting reproduction level
#' reproduction_level(NS_params)
#'
#' # The reproduction level can be changed without changing the steady state:
#' reproduction_level(params) <- 0.9
#' reproduction_level(params)
reproduction_level <- function(params) {
    UseMethod("reproduction_level")
}
#' @export
reproduction_level.MizerParams <- function(params) {
    params <- validParams(params)
    definition <- tuning_definition(params, "reproduction")
    if (is.null(definition) || is.null(definition$level_fun) ||
            !level_binding_current(params, definition)) {
        stop("Missing current law-specific reproduction-level support for ",
            "the effective reproduction provider.",
            call. = FALSE
        )
    }
    rates <- tuning_rate_context(params, params@initial_effort, 0)
    args <- tuning_context_args(params, rates, 0, params@initial_effort)
    value <- invoke_tuning_callback(
        definition$level_fun, params, args,
        "Reproduction-level callback"
    )
    if (!is.numeric(value) || is.complex(value) ||
            length(value) != nrow(params@species_params) ||
            any(is.nan(value)) || any(is.infinite(value))) {
        stop("The reproduction-level callback must return one numeric value ",
            "or NA for each species.",
            call. = FALSE
        )
    }
    species <- as.character(params@species_params$species)
    if (!is.null(names(value))) {
        if (anyNA(names(value)) || anyDuplicated(names(value)) ||
                !setequal(names(value), species)) {
            stop("The reproduction-level callback must name every species ",
                "exactly once.",
                call. = FALSE
            )
        }
        value <- value[species]
    }
    stats::setNames(as.numeric(value), species)
}
#' @rdname setBevertonHolt
#' @param value A number between 0 and 1, or a vector of numbers, giving the
#'   reproduction level for each species.
#' @export
`reproduction_level<-` <- function(params, value) {
    balanceRDD(params, reproduction_level = value)
}
#' Alias for `setBevertonHolt()`
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' An alias provided for backward compatibility with mizer version <= 2.0.4
#' @inherit setBevertonHolt
#' @export
#' @concept deprecated
setRmax <- setBevertonHolt
