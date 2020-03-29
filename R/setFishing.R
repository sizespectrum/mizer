#' Set fishing parameters
#' 
#' @section Setting fishing:
#' 
#' \strong{Gears}
#' 
#' In `mizer`, fishing mortality is imposed on species by fishing gears. The
#' total fishing mortality is obtained by summing over the mortality from all
#' gears,
#' \deqn{\mu_{f.i}(w) = \sum_g F_{g,i}(w),}
#' where the fishing mortality \eqn{F_{g,i}(w)} imposed by gear \eqn{g} on
#' species \eqn{i} at size \eqn{w} is calculated as:
#' \deqn{F_{g,i}(w) = S_{g,i}(w) Q_{g,i} E_{g},}
#' where \eqn{S} is the selectivity by species, gear and size, \eqn{Q} is the 
#' catchability by species and gear and \eqn{E} is the fishing effort by gear.
#' 
#' At the moment a species can only be selected by one fishing gear, although 
#' each gear can select more than one species (this is a limitation with the 
#' current package that will be developed in future releases). The gear
#' selecting each species can be specified in the `gear` column in the
#' species_params data frame. If no gear is specified, the default gear is
#' "knife_edge_gear".
#' 
#' \strong{Selectivity}
#' 
#' The selectivity at size of each gear has a range between 0 (not selected at
#' that size) to 1 (fully selected at that size). It is given by a selectivity
#' function. The name of the selectivity function is given by the `sel_func`
#' column in the species parameters data frame. Some selectivity functions are
#' included in the package: `knife_edge()`, `sigmoid_length()`,
#' `double_sigmoid_length()`, and `sigmoid_weight()`. New functions can be defined 
#' by the user. Each
#' gear has the same selectivity function for all the species it selects, but
#' the parameter values for each species may be different, e.g. the lengths of
#' species that a gear selects may be different.
#' 
#' Each selectivity function has a range of parameters. Values for these
#' parameters must be included as columns in the species parameters data.frame.
#' The names of the columns must exactly match the names of the corresponding
#' arguments of the selectivity function. For example, the default selectivity
#' function is `knife_edge()` which has sudden change of selectivity from 0 to 1
#' at a certain size. In its help page you can see that the `knife_edge()`
#' function has arguments `w` and `knife_edge_size` The first argument, `w`, is
#' size (the function calculates selectivity at size). All selectivity functions
#' must have `w` as the first argument. The values for the other arguments must
#' be found in the species parameters data.frame. So for the `knife_edge()`
#' function there should be a `knife_edge_size` column. Because `knife_edge()`
#' is the default selectivity function, the `knife_edge_size` argument has a
#' default value = `w_mat`.
#' 
#' \strong{Catchability}
#' 
#' Catchability is used as an additional scalar to make the link between gear 
#' selectivity, fishing effort and fishing mortality. For example, it can be set 
#' so that an effort of 1 gives a desired fishing mortality.
#' In this way effort can then be specified relative to a 'base effort', e.g. 
#' the effort in a particular year.
#' 
#' Because currently mizer only allows one gear to select each species, the
#' catchability matrix \eqn{Q_{g,i}} reduces to a catchability vector 
#' \eqn{Q_{i}}, and this is given as a column `catchability` in the species
#' parameter data frame. If it is not specified, it defaults to 1.
#' 
#' \strong{Effort}
#' 
#' The initial fishing effort is stored in the `MizerParams` object. If it is
#' not supplied, it is set to zero. The initial effort can be overruled when
#' the simulation is run with `project()`, where it is also possible to specify
#' an effort that varies through time.
#' 
#' @param params A MizerParams object
#' @param selectivity An array (gear x species x w) that holds the selectivity of
#'   each gear for species and size, \eqn{S_{g,i,w}}.
#' @param catchability An array (gear x species) that holds the catchability of
#'   each species by each gear, \eqn{Q_{g,i}}.
#' @param initial_effort Optional. A number or a named numeric vector specifying
#'   the fishing effort. If a number, the same effort is used for all gears. If
#'   a vector, must be named by gear.
#' @param ... Unused
#'   
#' @return MizerParams object with updated catchability and selectivity. Because
#'   of the way the R language works, `setFishing()` does not make the changes
#'   to the params object that you pass to it but instead returns a new params
#'   object. So to affect the change you call the function in the form
#'   `params <- setFishing(params, ...)`.
#' @export
#' @family functions for setting parameters
#' @examples
#' \dontrun{
#' params <- NS_params
#' # Change knife edge size for species 1
#' params@species_params$knife_edge_size[1] <- 15
#' params <- setFishing(params)
#' }
setFishing <- function(params, selectivity = NULL, catchability = NULL, 
                       initial_effort = NULL, ...) {
    assert_that(is(params, "MizerParams"))
    species_params <- params@species_params
    gear_params <- params@gear_params
    no_sp <- nrow(species_params)
    no_gears <- length(unique(gear_params$gear))
    no_w <- length(params@w)
    
    if (!is.null(selectivity)) {
        assert_that(length(dim(selectivity)) == 3,
                    dim(selectivity)[[2]] == no_sp,
                    dim(selectivity)[[3]] == length(params@w))
        params@selectivity <- selectivity
    } else {
        selectivity <- 
            array(0, dim = c(no_gears, no_sp, no_w),
                  dimnames = list(gear = as.character(unique(gear_params$gear)), 
                                  sp = as.character(params@species_params$species),
                                  w = signif(params@w, 3)
                  )
            )
        for (g in seq_len(nrow(gear_params))) {
            # get args
            # These as.characters are annoying - but factors everywhere
            arg <- names(formals(as.character(gear_params[g, 'sel_func'])))
            # lop off w as that is always the first argument of the selectivity functions
            arg <- arg[!(arg %in% "w")]
            if (!all(arg %in% colnames(gear_params))) {
                stop("Some arguments needed for the selectivity function are ",
                     "missing in the gear_params dataframe.")
            }
            # Check that there are no missing values for selectivity parameters
            if (any(is.na(as.list(gear_params[g, arg])))) {
                stop("Some selectivity parameters are NA.")
            }
            # Call selectivity function with selectivity parameters
            par <- c(w = list(params@w), as.list(gear_params[g, arg]))
            sel <- do.call(as.character(gear_params[g, 'sel_func']), args = par)
            if (!is.null(comment(params@selectivity)) &&
                any(params@selectivity[as.character(species_params[g, 'gear']),
                                       as.character(gear_params[g, 'species']),
                ] != sel)) {
                message("The selectivity has been commented and therefore will ",
                        "not be recalculated from the species parameters.")
            } else {
                params@selectivity[as.character(gear_params[g, 'gear']), 
                                   as.character(gear_params[g, 'species']),
                ] <- sel
            }
        }
    }
    
    if (!is.null(catchability)) {
        assert_that(length(dim(catchability)) == 2,
                    dim(selectivity)[[2]] == no_sp)
        params@catchability <- catchability
    } else {
        # If no catchability column in species_params, set to 1
        species_params <- set_species_param_default(species_params,
                                                    "catchability", 1)
        
        for (g in seq_len(no_sp)) {
            # Now do catchability
            if (!is.null(comment(params@catchability)) &&
                any(params@catchability[as.character(species_params[g,'gear']), g] !=
                    species_params[g, "catchability"])) {
                message("The catchability has been commented and therefore will ",
                        "not be updated from the species parameters.")
            } else {
                params@catchability[as.character(species_params[g,'gear']), g] <- 
                    species_params[g, "catchability"]
            }
        }
    }
    
    if (!is.null(initial_effort)) {
        validate_effort_vector(params, initial_effort)
        params@initial_effort[] <- initial_effort
        comment(params@initial_effort) <- comment(initial_effort)
    }
    
    params@species_params <- species_params
    return(params)
}

#' @rdname setFishing
#' @export
gear_params <- function(params) {
    params@gear_params
}

#' @rdname setFishing
#' @param value A data frame with the species parameters
#' @export
`gear_params<-` <- function(params, value) {
    value <- validGearParams(value, params@species_params)
    params@gear_params <- value
    setFishing(params)
}

#' @rdname setFishing
#' @export
getCatchability <- function(params) {
    params@catchability
}

#' @rdname setFishing
#' @export
getSelectivity <- function(params) {
    params@selectivity
}

#' @rdname setFishing
#' @export
getInitialEffort <- function(params) {
    params@initial_effort
}

#' Check validity of gear parameters and set defaults for missing but
#' required parameters or transfer them from species_params if available
#' 
#' @param gear_params Gear parameter data frame
#' @param species_params Species parameter data frame
#' @return A valid gear parameter data frame
validGearParams <- function(gear_params, species_params) {
    assert_that(is.data.frame(gear_params),
                is.data.frame(species_params))
    
    no_sp <- nrow(species_params)
    
    if (nrow(gear_params) < 1) {
        if (!is.null(species_params$gear)) {
            # Try to take parameters from species_params
            gear_params <- 
                data.frame(gear = species_params$gear,
                           species = species_params$species)
            if (!is.null(species_params$sel_func)) {
                gear_params$sel_func <- species_params$sel_func
            } else {
                gear_params$sel_func <- "knife_edge"
            }
            # copy over any selectivity function parameters
            for (g in seq_len(no_sp)) {
                args <- names(formals(as.character(gear_params[g, 'sel_func'])))
                args <- args[!(args %in% "w")]
                for (arg in args) {
                    if (!is.null(species_params[[arg]])) {
                        gear_params[[arg]] <- species_params[[arg]]
                    } else if (arg == "knife_edge_size") {
                        gear_params[[arg]] = species_params$w_mat
                    }
                }
            }
        } else {
            gear_params <- 
                data.frame(species = species_params$species,
                           gear = "knife_edge_gear",
                           sel_func = "knife_edge",
                           knife_edge_size = species_params$w_mat)
        }
    }
    
    # TODO: Check that there are no duplicate gear-species pairs
    
    # Check that every row is complete
    for (g in 1:nrow(gear_params)) {
        # get args
        # These as.characters are annoying - but factors everywhere
        arg <- names(formals(as.character(gear_params[g, 'sel_func'])))
        # lop off w as that is always the first argument of the selectivity functions
        arg <- arg[!(arg %in% "w")]
        if (!all(arg %in% colnames(gear_params))) {
            stop("Some arguments needed for the selectivity function are ",
                 "missing in the gear parameter dataframe.")
        }
        # Check that there are no missing values for selectivity parameters
        if (any(is.na(as.list(gear_params[g, arg])))) {
            stop("Some selectivity parameters are NA.")
        }
    }
    gear_params
}

#' Check that an effort vector is specified correctly
#' 
#' Throws an error with an explanatory message when the supplied `effort`
#' vector is not valid for the model described by `params`.
#' 
#' @param params A MizerParams object
#' @param effort An effort vector
#' 
#' @return TRUE if `effort` is valid. Throws an error otherwise.
#' @export
#' @keywords internal
#' @concept helper
validate_effort_vector <- function(params, effort) {
    assert_that(is(params, "MizerParams"),
                is.numeric(effort))
    no_gears <- dim(params@catchability)[1]
    if ((length(effort) > 1) & (length(effort) != no_gears)) {
        stop("Effort vector must be the same length as the number of fishing gears\n")
    }
    # If more than 1 gear need to check that gear names match
    gear_names <- dimnames(params@catchability)[[1]]
    effort_gear_names <- names(effort)
    if (length(effort) == 1 & is.null(effort_gear_names)) {
        effort_gear_names <- gear_names
    }
    if (!all(gear_names %in% effort_gear_names)) {
        stop("Gear names in the MizerParams object (", 
             paste(gear_names, collapse = ", "), 
             ") do not match those in the effort vector.")
    }
    return(TRUE)
}