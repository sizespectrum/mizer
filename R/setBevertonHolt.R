#' Set Beverton-Holt density dependence
#'
#' `r lifecycle::badge("experimental")`
#' Takes a MizerParams object `params` with arbitrary density dependence and
#' returns a MizerParams object with Beverton-Holt density-dependence in such a
#' way that the energy invested into reproduction by the mature individuals
#' leads to the reproduction rate that is required to maintain the given egg
#' abundance. Hence if you have tuned your `params` object to describe a
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
#' \deqn{R_{di} = \frac{\rm{erepro}}{2 w_{min}} E_R,}{R_di = (erepro/(2 w_min)) E_R,}
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
#' ```{r Beverton-Holt-plot, echo = FALSE, fig.height = 2.5, fig.width = 5}
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
#'                   colour = dd, size = dd)) +
#'     geom_point(aes(x = 5/4, y = 5/6), size = 2) +
#'     labs(linetype = "", size = "R_max", colour = "R_max") +
#'     scale_y_continuous(name = "Reproduction rate [eggs/year]",
#'                        breaks = c(5/6), labels = c("R_dd")) +
#'     scale_x_continuous(name = "Energy invested [g/year]",
#'                        breaks = c(5/4), labels = c("E_R")) +
#'     scale_size_manual(values = c(1, 0.5)) +
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
#' calculate the energy \eqn{E_R} that is invested into reproduction by the
#' mature individuals and the reproduction rate \eqn{R_{dd}}{R_dd} that is
#' required to keep the egg abundance constant. These two values determine the
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
#' up to `Inf`. If a smaller value is requested a warning is issued and the
#' value is increased to the value required for a reproduction level of 0.99.
#' 
#' The values for the `reproduction_level` must be positive and
#' less than 1. The values for `erepro` must be large enough to allow the
#' required reproduction rate. If a smaller value is requested a warning is
#' issued and the value is increased to the smallest possible value. The values
#' for `erepro` should also be smaller than 1 to be physiologically sensible,
#' but this is not enforced by the function.
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
#' @param R_factor `r lifecycle::badge("deprecated")` Use
#'   `reproduction_level = 1 / R_factor` instead.
#'
#' @return A MizerParams object
#' @examples
#' params <- NS_params
#' species_params(params)$erepro
#' # Attempting to set the same erepro for all species
#' params <- setBevertonHolt(params, erepro = 0.1)
#' t(species_params(params)[, c("erepro", "R_max")])
#' # Setting erepro for some species
#' params <- setBevertonHolt(params, erepro = c("Gurnard" = 0.6, "Plaice" = 0.95))
#' t(species_params(params)[, c("erepro", "R_max")])
#' # Setting R_max
#' R_max <- 1e17 * species_params(params)$w_inf^-1
#' params <- setBevertonHolt(NS_params, R_max = R_max)
#' t(species_params(params)[, c("erepro", "R_max")])
#' # Setting reproduction_level
#' params <- setBevertonHolt(params, reproduction_level = 0.3)
#' t(species_params(params)[, c("erepro", "R_max")])
#' @export
setBevertonHolt <- function(params, R_factor = deprecated(), erepro, 
                            R_max, reproduction_level) {
    assert_that(is(params, "MizerParams"))
    no_sp <- nrow(params@species_params)
    
    # check number of arguments
    num_args <- hasArg("erepro") + hasArg("R_max") + 
        hasArg("reproduction_level") + hasArg("R_factor")
    if (num_args > 1) {
        stop("You should only provide `params` and one other argument.")
    }
    if (num_args == 0) {
        # no values given, so use previous erepro
        erepro <- species_params(params)$erepro
    }
    
    # No matter which argument is given, I want to manipulate the values
    if (!missing("erepro")) values <- erepro
    if (hasArg("R_max")) values <- R_max
    if (hasArg("reproduction_level")) values <- reproduction_level
    if (hasArg("R_factor")) values <- R_factor
    
    if (length(values) == 1 && is.null(names(values))) {
        values <- rep(values, no_sp)
    }
    if (length(values) != no_sp && is.null(names(values))) {
        stop("You need to supply a vector of length ", no_sp,
             " or a single number or a named vector.")
    }
    if (is.null(names(values))) {
        names(values) <- params@species_params$species
    }
    values <- values[!is.na(values)]
    if (length(values) == 0) return(params) # Nothing to do
    if (!all(is.numeric(values))) {
        stop("You provided invalid non-numeric values.")
    }
    # select the species that are affected
    species <- valid_species_arg(params, names(values))
    sp_idx <- match(species, params@species_params$species)

    rdi <- getRDI(params)[species]
    if (any(rdi == 0)) { # This should never happen, but did happen in the past.
        stop("Some species have no reproduction.")
    }
    params@rates_funcs$RDD <- "BevertonHoltRDD"
    
    rdd_new <- getRequiredRDD(params)[species]
    
    if (!missing(erepro)) {
        erepro_new <- values
        erepro_old <- params@species_params$erepro[sp_idx]
        rdi_new <- rdi * erepro_new / erepro_old
        wrong <- rdi_new < rdd_new
        if (any(wrong)) {
            rdi_new[wrong] <- rdd_new[wrong]
            erepro_new[wrong] <- (erepro_old * rdi_new / rdi)[wrong]
            warning("For the following species `erepro` ",
                    "has been increased to the smallest ",
                    "possible value: ",
                    paste0("erepro[", species[wrong], "] = ",
                           signif(erepro_new[wrong], 3),
                           collapse = "; "))
        }
        r_max_new <- rdi_new * rdd_new / (rdi_new - rdd_new)
        r_max_new[is.nan(r_max_new)] <- Inf
        params@species_params$erepro[sp_idx] <- erepro_new
        params@species_params$R_max[sp_idx] <- r_max_new
        
        params@time_modified <- lubridate::now()
        return(params)
    }
    
    # We now know that we are setting R_max, which however can have been
    # specified in different ways.
    if (!missing(reproduction_level)) {
        if (!all(values >= 0 & values < 1)) {
            stop("The reproduction level must be smaller than 1 and non-negative.")
        }
        r_max_new <- rdd_new / values
    }
    if (!missing(R_factor)) {
        if (!all(values > 1)) {
            stop("The R_factor must be greater than 1.")
        }
        r_max_new <- rdd_new * values
    }
    if (!missing(R_max)) {
        wrong <- values < rdd_new
        if (any(wrong)) {
            warning("For the following species the requested `R_max` ",
                    "was too small and has been increased to give a ",
                    "reproduction level of 0.99: ",
                    paste(species[wrong], collapse = ", "))
            values[wrong] <- rdd_new[wrong] / 0.99
        }
        r_max_new <- values
    }
    
    # Be careful to choose the expression that works also with r_max_new = Inf
    # rdi_new <- r_max_new * rdd_new / (r_max_new - rdd_new)
    rdi_new <- rdd_new / (1 - rdd_new / r_max_new)
    
    params@species_params$R_max[sp_idx] <- r_max_new
    params@species_params$erepro[sp_idx] <-
        params@species_params$erepro[sp_idx] * rdi_new / rdi
    wrong <- params@species_params$erepro[sp_idx] > 1
    if (any(wrong)) {
        warning("The following species require an unrealistic reproductive ",
                "efficiency greater than 1: ",
                paste(species[wrong], collapse = ", "))
    }
    
    params@time_modified <- lubridate::now()
    return(params)
}

getRequiredRDD <- function(params) {
    # Calculate required rdd
    mumu <- getMort(params)
    gg <- getEGrowth(params)
    rdd_new <- getRDD(params) # to get the right structure
    for (i in seq_len(nrow(params@species_params))) {
        gg0 <- gg[i, params@w_min_idx[i]]
        if (!(gg0 > 0)) {
            warning("Eggs of species ", params@species_params$species[i],
                    " have zero growth rate.")
        }
        mumu0 <- mumu[i, params@w_min_idx[i]]
        DW <- params@dw[params@w_min_idx[i]]
        n0 <- params@initial_n[i, params@w_min_idx[i]]
        if (!(n0 > 0)) {
            warning("Species ", params@species_params$species[i],
                    "appears to have no eggs.")
        }
        rdd_new[i] <- n0 * (gg0 + DW * mumu0)
    }
    rdd_new
}

#' Get reproduction level
#'
#' `r lifecycle::badge("experimental")`
#' The reproduction level is the ratio between the density-dependent
#' reproduction rate and the maximal reproduction rate.
#'
#' @param params A MizerParams object
#'
#' @return A named vector with the reproduction level for each species.
#' @export
#' @examples
#' getReproductionLevel(NS_params)
#'
#' # The reproduction level can be changed without changing the steady state:
#' params <- setBevertonHolt(NS_params, reproduction_level = 0.9)
#' getReproductionLevel(params)
#'
#' # The result is the ratio of RDD and R_max
#' identical(getRDD(params) / species_params(params)$R_max,
#'           getReproductionLevel(params))
getReproductionLevel <- function(params) {
    assert_that(is(params, "MizerParams"))
    if (!"R_max" %in% names(params@species_params)) {
        stop("No `R_max` is included in the species parameters.")
    }
    getRDD(params) / params@species_params$R_max
}


#' Alias for `setBevertonHolt()`
#' 
#' @description
#' `r lifecycle::badge("deprecated")`
#' 
#' An alias provided for backward compatibility with mizer version <= 2.0.4
#' @inherit setBevertonHolt
#' @export
#' @concept deprecated
setRmax <- setBevertonHolt
