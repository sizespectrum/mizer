#' Check whether two objects are different
#'
#' Check whether two objects are numerically different, ignoring all attributes.
#'
#' We use this helper function in particular to see if a new value for a slot
#' in MizerParams is different from the existing value in order to give the
#' appropriate messages.
#'
#' @param a First object
#' @param b Second object
#'
#' @return TRUE or FALSE
#' @concept helper
different <- function(a, b) {
    !isTRUE(all.equal(a, b, check.attributes = FALSE, scale = 1,
                      tolerance = 10 * .Machine$double.eps))
}

#' Extract a column from a parameter table without partial matching
#'
#' The `$` operator on a data frame partially matches column names, so
#' `species_params(params)$a` would return the `alpha` column when there is no
#' `a` column. That silently gives the wrong parameter, so the `$` methods for
#' the mizer parameter tables use this helper instead, which matches exactly.
#'
#' When the name would have partially matched a single column under the old
#' behaviour we warn, because that is precisely the case where existing user
#' code changes its meaning. A name that matches nothing at all returns `NULL`
#' silently, so that `is.null(species_params(params)$foo)` stays a usable way
#' of testing whether a parameter is present.
#'
#' @param x The parameter table.
#' @param name The name of the column.
#' @param what The name of the class, used in the warning message.
#'
#' @return The column, or `NULL` if there is no column of that name.
#' @concept helper
#' @noRd
exact_column <- function(x, name, what) {
    out <- .subset2(x, name, exact = TRUE)
    if (is.null(out)) {
        partial <- names(x)[startsWith(names(x), name)]
        if (length(partial) == 1) {
            warning("There is no `", name, "` column in this `", what,
                    "` table, so `$", name, "` returns NULL. Earlier versions ",
                    "of mizer partially matched the column name and returned ",
                    "the `", partial, "` column instead. Use `$", partial,
                    "` if that is what you wanted.")
        }
    }
    out
}

#' Bin average of a power law over geometric bins
#'
#' Computes the exact average of the power law \eqn{w^d} over each bin
#' \eqn{[w_j, w_{j+1}]}, i.e.
#' \deqn{\overline{w^d}_j = \frac{1}{\Delta w_j}\int_{w_j}^{w_{j+1}} w^d\, dw.}
#'
#' The integral has a closed form, so the result is exact (not merely second
#' order):
#' \deqn{\overline{w^d}_j = \frac{w_{j+1}^{d+1} - w_j^{d+1}}{(d+1)\,\Delta w_j},
#'   \quad d \neq -1,}
#' \deqn{\overline{w^d}_j = \frac{\ln(w_{j+1}/w_j)}{\Delta w_j},
#'   \quad d = -1.}
#'
#' This is used by the bin-averaged (second-order) code paths that need the
#' average of a power-law rate over each bin, for example [setExtMort()] and
#' the resource capacity and rate in [setResource()]. The grid does not need to
#' be geometric; only the left bin edges `w` and the bin widths `dw` are used,
#' with \eqn{w_{j+1} = w_j + \Delta w_j}.
#'
#' An optional upper cutoff `w_max` handles a knife-edge truncation of the power
#' law (for example the resource carrying capacity, which is \eqn{\kappa
#' w^{-\lambda}} below `w_pp_cutoff` and zero above it). The bin straddling the
#' cutoff then receives the *partial* bin-average — the power-law average over
#' the part of the bin below `w_max`, divided by the full bin width — and bins
#' entirely above the cutoff get zero.
#'
#' @param w Numeric vector of left bin edges \eqn{w_j}.
#' @param dw Numeric vector of bin widths \eqn{\Delta w_j} (same length as `w`).
#' @param d Single numeric exponent of the power law.
#' @param w_max Optional upper cutoff. The power law is taken to be zero above
#'   `w_max`, with the straddling bin getting the partial average. Defaults to
#'   `Inf` (no cutoff), in which case the result is identical to the uncut
#'   formula.
#'
#' @return A numeric vector (same length as `w`) of bin averages of
#'   \eqn{w^d} (truncated at `w_max` when supplied).
#' @concept helper
#' @keywords internal
power_law_bin_average <- function(w, dw, d, w_max = Inf) {
    # Upper integration edge of each bin, clipped at the cutoff. The pmax keeps
    # bins lying entirely above the cutoff from contributing a negative (they
    # collapse to a zero-width interval and so average to zero).
    w_next <- pmax(pmin(w + dw, w_max), w)
    if (d == -1) {
        log(w_next / w) / dw
    } else {
        (w_next^(d + 1) - w^(d + 1)) / ((d + 1) * dw)
    }
}

#' Construct the background resource power-law spectrum
#'
#' Internal helper returning the auto-calculated resource power law
#' \eqn{\kappa\, w^{-\lambda}} on the full size grid, optionally truncated at an
#' upper cutoff `w_max`. When the `bin_average` entry of the model's
#' `second_order_w` slot is set, the exact bin average of the power law over
#' each bin is returned instead (with the bin straddling `w_max` getting the
#' partial average), so that the initial resource is the finite-volume cell
#' average of the background spectrum, consistent with the bin-integrated
#' encounter convolution that consumes it as a cell average. Otherwise the
#' left-edge point values are returned, byte-identical to previous mizer.
#'
#' This is used wherever the background resource spectrum \eqn{\kappa
#' w^{-\lambda}} is constructed from scratch: the initial resource abundance in
#' [newMultispeciesParams()] and the temporary prey spectra used to compute the
#' default `gamma`/`f0` ([get_gamma_default()], [get_f0_default()]) and the
#' consumer initial abundances ([get_initial_n()]). The bin-averaged resource
#' capacity and rate are produced directly by [setResource()].
#'
#' @param params A MizerParams object whose `second_order_w` slot controls the
#'   gating and whose `w_full`/`dw_full` give the grid.
#' @param kappa The coefficient \eqn{\kappa} of the power law.
#' @param lambda The exponent so the power law is \eqn{w^{-\lambda}}.
#' @param w_max Optional upper cutoff. The power law is taken to be zero at and
#'   above `w_max`. Defaults to `Inf` (no cutoff).
#'
#' @return A numeric vector (same length as `w_full`) of the resource spectrum.
#' @concept helper
#' @keywords internal
resource_power_law <- function(params, kappa, lambda, w_max = Inf) {
    w_full <- params@w_full
    if (isTRUE(params@second_order_w[["bin_average"]])) {
        kappa * power_law_bin_average(w_full, params@dw_full, -lambda,
                                      w_max = w_max)
    } else {
        n_pp <- kappa * w_full ^ (-lambda)
        n_pp[w_full >= w_max] <- 0
        n_pp
    }
}

#' Geometric bin centres of the size grid
#'
#' Internal helper for the second-order plotting code. A finite-volume cell
#' average \eqn{N_j = (1/\Delta w_j)\int_{w_j}^{w_{j+1}} N\,dw} does not live at
#' the left bin edge \eqn{w_j} but at the geometric bin centre
#' \deqn{w^*_j = \sqrt{w_j\,w_{j+1}} = w_j\sqrt{\beta},}
#' where \eqn{\beta = w_{j+1}/w_j} is the (constant) bin ratio of the
#' logarithmic grid. This is the log-symmetric, second-order-correct location at
#' which to plot a bin-averaged quantity (it is exact for the community spectrum
#' \eqn{N\propto w^{-2}}). It is a uniform half-bin shift to the right on the log
#' axis, the same for the consumer grid `w` and the full prey/resource grid
#' `w_full`.
#'
#' @param params A MizerParams object.
#' @param w_full If `TRUE`, return the centres of the full (resource) grid
#'   `params@w_full`; otherwise the consumer grid `params@w`.
#' @return A numeric vector of geometric bin centres, one per grid node.
#' @concept helper
#' @keywords internal
bin_midpoints <- function(params, w_full = FALSE) {
    # The grid is geometric, so the bin ratio beta is constant across the whole
    # (consumer and resource) grid.
    beta <- params@w_full[2] / params@w_full[1]
    w <- if (w_full) params@w_full else params@w
    w * sqrt(beta)
}

#' Length-weight conversion
#'
#' For each species, convert between length and weight using the relationship
#' \deqn{w_i = a_i l_i^{b_i}}{w_i = a_i l_i^b_i} or
#' \deqn{l_i = (w_i / a_i)^{1/b_i}}{l_i = (w_i / a_i)^{1/b_i}}
#' where `a` and `b` are taken from the species parameter data frame and
#' \eqn{i}{i} is the species index.
#'
#' This is useful for converting a length-based species parameter to a
#' weight-based species parameter.
#'
#' @param l Lengths in cm. Either a single number used for all species or a
#'   vector with one number for each species.
#' @param w Weights in grams. Either a single number used for all species or a
#'   vector with one number for each species.
#' @param species_params A species parameter data frame or a MizerParams object.
#' @return A vector with one entry for each species. `l2w()` returns a vector
#' of weights in grams and `w2l()` returns a vector of lengths in cm.
#' @export
#' @concept helper
l2w <- function(l, species_params) {
    # The checks avoid `assert_that()`, which builds its message whether or not
    # it is needed.
    if (!is.numeric(l)) stop("l is not a numeric or integer vector")
    sp <- species_params
    if (inherits(species_params, "MizerParams")) {
        sp <- species_params@species_params
    }
    if (!is.data.frame(sp)) {
        stop("The second argument must be either a MizerParams object or a
             species paramter data frame.")
    }
    if (length(l) != 1 && length(l) != nrow(sp)) {
        stop("The length of 'l' should be one or equal to the number of species.")
    }
    if (!all(c("a", "b") %in% names(sp))) {
        stop("The species parameter data frame must contain columns 'a' and 'b'.")
    }

    sp[["a"]] * l^sp[["b"]]
}

#' @rdname l2w
#' @export
w2l <- function(w, species_params) {
    if (!is.numeric(w)) stop("w is not a numeric or integer vector")
    sp <- species_params
    if (inherits(species_params, "MizerParams")) {
        sp <- species_params@species_params
    }
    if (!is.data.frame(sp)) {
        stop("The second argument must be either a MizerParams object or a
             species paramter data frame.")
    }
    if (length(w) != 1 && length(w) != nrow(sp)) {
        stop("The length of 'w' should be one or equal to the number of species.")
    }
    if (!all(c("a", "b") %in% names(sp))) {
        stop("The species parameter data frame must contain columns 'a' and 'b'.")
    }

    (w / sp[["a"]])^(1 / sp[["b"]])
}

#' Calculate steady state abundance
#'
#' Solves the steady-state transport equation with growth, mortality and
#' diffusion held fixed and the supplied abundance in each egg size class.
#'
#' @param params A MizerParams object
#' @param g A matrix of growth rates (species x size)
#' @param mu A matrix of mortality rates (species x size)
#' @param D A matrix of diffusion rates (species x size)
#' @param N0 A vector with the abundance at the smallest size for each species
#' @param max_iterations Maximum number of defect-correction updates for a
#'   van Leer flux limiter.
#' @param tol Relative tolerance for both the correction and the transport
#'   residual, checked separately for every species.
#' @param relax Under-relaxation factor in (0, 1] for defect correction.
#' @return A matrix with the steady-state abundance. Stops with an error naming
#'   the unconverged species if the iteration limit is reached.
#'
#' @details
#' The flux scheme is read from [second_order_w()]. The upwind and centred
#' systems are linear and are solved directly. With a van Leer limiter, defect
#' correction solves the nonlinear equation using an upwind preconditioner with
#' the same face placement and diffusion as the second-order scheme. Only the
#' preconditioner has its reconstruction weights set to zero: the residual
#' retains the model's chosen flux limiter.
#'
#' Each present species is solved with unit egg abundance and then rescaled to
#' `N0`, because the limiter is unchanged by multiplying a species' spectrum by
#' a positive constant. This prevents species with large egg abundances from
#' hiding non-convergence in other species. Convergence requires both a small
#' correction and a small residual of the original transport equation, so
#' clipping negative iterates cannot turn a stalled iteration into a solution.
#'
#' Abundance is held at zero below the egg size and above the support top given
#' by `support_top_idx()`, matching the boundary conditions of [project()].
#' @concept helper
get_steady_state_n <- function(params, g, mu, D, N0,
                               max_iterations = 500, tol = 1e-10,
                               relax = 0.5) {
    flux_limiter <- flux_limiter_scheme(params)
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    zero <- matrix(0, nrow = no_sp, ncol = no_w,
                   dimnames = list(params@species_params$species,
                                   dimnames(params@initial_n)[[2]]))
    j_start <- params@w_min_idx
    idxs_start <- cbind(seq_len(no_sp), j_start)
    w_top <- support_top_idx(params)

    # Set the egg density rather than a recruitment flux. The time-stepping
    # diagonal contains the identity, which is absent from the steady equation.
    steady_coefs <- function(coefs) {
        coefs$b <- coefs$b - 1
        coefs$a[idxs_start] <- 0
        coefs$b[idxs_start] <- 1
        coefs$c[idxs_start] <- 0
        coefs
    }
    solve_coefs <- function(coefs, rhs) {
        project_n_loop(zero, coefs$a, coefs$b, coefs$c, rhs, j_start)
    }

    if (flux_limiter == "none") {
        # Preserve the legacy solve, including its unnormalised right-hand
        # side, so existing first-order models keep exactly the same results.
        coefs <- get_transport_coefs(params, zero, g, mu, dt = 1,
                                     recruitment_flux = rep(0, no_sp),
                                     d = D, flux_limiter = flux_limiter)
        coefs <- steady_coefs(coefs)
        rhs <- zero
        rhs[idxs_start] <- N0
        return(zero_above_support(solve_coefs(coefs, rhs), w_top))
    }

    if (all(N0 == 0)) {
        return(zero)
    }
    present <- N0 > 0
    active <- col(zero) >= j_start & col(zero) <= w_top & present
    # Inactive cells have the equation n = 0, also preventing a singular or
    # non-finite inactive tail from contaminating the tridiagonal solve.
    active_coefs <- function(coefs) {
        coefs <- steady_coefs(coefs)
        coefs$a[!active] <- 0
        coefs$b[!active] <- 1
        coefs$c[!active] <- 0
        coefs
    }
    fail <- function(species) {
        species <- params@species_params$species[species]
        text <- paste0("The steady-state calculation did not converge for: ",
                       toString(species), ".")
        condition <- list(message = text, call = NULL, species = species)
        class(condition) <- c("mizer_steady_nonconvergence", "error",
                              "condition")
        stop(condition)
    }
    row_max <- function(x) apply(x, 1, max)
    check_finite_rows <- function(x) {
        bad <- rowSums(!is.finite(x)) > 0
        if (any(bad)) {
            fail(which(bad))
        }
    }

    rhs <- zero
    rhs[idxs_start] <- as.numeric(present)
    # Centred reconstruction is linear, so solve it directly. For van Leer,
    # do not use get_transport_coefs(..., flux_limiter = "none") here: that
    # selects the legacy scheme, whose fluxes sit at different faces. The
    # log-FV operator with chi = 0 is the appropriate upwind preconditioner.
    preconditioner <- if (flux_limiter == "centred") "centred" else "none"
    low <- get_transport_coefs_logfv(params, zero, g, mu, dt = 1,
                                     recruitment_flux = rep(0, no_sp),
                                     d = D, flux_limiter = preconditioner)
    low <- active_coefs(apply_upper_cutoff(low, w_top))
    n <- pmax(solve_coefs(low, rhs), 0)

    # Check the initial solution and then at most max_iterations updates.
    for (iteration in seq_len(max_iterations + 1L)) {
        # Check before evaluating chi: non-finite differences would otherwise
        # cause an unrelated NA-subscript error in the limiter.
        check_finite_rows(n)
        coefs <- get_transport_coefs(params, n, g, mu, dt = 1,
                                     recruitment_flux = rep(0, no_sp),
                                     d = D, flux_limiter = flux_limiter)
        coefs <- active_coefs(coefs)
        left <- coefs$a * cbind(0, n[, -no_w, drop = FALSE])
        centre <- coefs$b * n
        right <- coefs$c * cbind(n[, -1, drop = FALSE], 0)
        residual <- left + centre + right - rhs
        check_finite_rows(residual)
        correction <- solve_coefs(low, residual)
        check_finite_rows(correction)

        change <- row_max(abs(correction)) / pmax(row_max(abs(n)), 1)
        scale <- abs(left) + abs(centre) + abs(right) + abs(rhs)
        # Scale each cell by the terms in its balance. A round-off floor per
        # species avoids demanding relative accuracy in a numerically empty
        # tail.
        floor <- .Machine$double.eps * pmax(row_max(scale), 1)
        error <- row_max(abs(residual) / pmax(scale, floor))
        converged <- change < tol & error < tol
        if (all(converged)) {
            return(n * N0)
        }
        if (flux_limiter == "centred") {
            # A clipped linear solution must also satisfy the original
            # equation; repeating the same solve cannot repair a negative root.
            break
        }
        if (iteration <= max_iterations) {
            correction[converged, ] <- 0
            n <- pmax(n - relax * correction, 0)
            n[idxs_start] <- rhs[idxs_start]
        }
    }
    fail(which(!converged))
}

#' Warn about column names that look like mis-spelled parameter names
#'
#' Compares the actual column names of a parameter data frame against the set of
#' recognised names and raises a single [signal_info()] report about any names
#' that are a near match (a likely typo) but not an exact match. Names are not
#' corrected, only flagged, so that legitimate custom columns are left
#' untouched.
#'
#' @param actual Character vector of the column names present in the data frame.
#' @param known Character vector of recognised parameter names.
#' @param df_type A short description of the data frame used in the message,
#'   e.g. `"gear parameter"` or `"species parameter"`.
#' @param var A string naming the quantity the report is about, see
#'   [signal_info()]. Given explicitly rather than derived from `df_type`, which
#'   is prose and may be reworded: `var` is the key that handlers and
#'   `with_info_level(except = )` match on.
#' @param curated Character vector of known problematic names that should always
#'   be flagged even when they are further than the fuzzy-match threshold from
#'   any recognised name (e.g. familiar abbreviations).
#'
#' @return Called for its side effect (an information signal). Returns `NULL`
#'   invisibly.
#' @concept helper
#' @noRd
check_for_misspellings <- function(actual, known, df_type, var,
                                   curated = character(0)) {
    unknown <- setdiff(actual, known)
    if (length(unknown) == 0) {
        return(invisible(NULL))
    }
    # Only fuzzy-match against names long enough that a near miss is meaningful.
    known_fuzzy <- known[nchar(known) >= 3]
    hits <- character(0)
    guess <- character(0)
    for (nm in unknown) {
        if (length(known_fuzzy) == 0) {
            j <- integer(0)
            near <- FALSE
        } else {
            d <- as.vector(adist(nm, known_fuzzy, ignore.case = TRUE))
            j <- which.min(d)
            near <- nchar(nm) >= 3 &&
                (d[j] == 1 || (d[j] <= 2 && nchar(known_fuzzy[j]) >= 6))
        }
        if (nm %in% curated || near) {
            hits <- c(hits, nm)
            guess <- c(guess, if (length(j) == 1) known_fuzzy[j] else NA_character_)
        }
    }
    if (length(hits) > 0) {
        suggestions <- ifelse(is.na(guess),
                              paste0("`", hits, "`"),
                              paste0("`", hits, "` (did you mean `", guess,
                                     "`?)"))
        signal_info(var,
                    paste0("Some column names in your ", df_type, " data ",
                           "frame are very close to standard parameter names: ",
                           paste(suggestions, collapse = ", "),
                           ". Please check for mis-spellings."),
                    level = 1, severity = "warning", unhandled = "show")
    }
    invisible(NULL)
}

#' Convert an S4 MizerParams or MizerSim object to an S3 object
#'
#' @param object An S4 MizerParams or MizerSim object.
#' @return An S3 MizerParams or MizerSim object (a list).
#' @keywords internal
upgrade_s4_to_s3 <- function(object) {
    if (!isS4(object)) return(object)
    attrs <- attributes(object)
    sl_names <- setdiff(names(attrs), c("class", "comment", "dotData"))
    obj_list <- lapply(sl_names, function(s) attrs[[s]])
    names(obj_list) <- sl_names
    cl <- class(object)
    attr(cl, "package") <- NULL
    class(obj_list) <- cl
    if ("comment" %in% names(attrs)) {
        comment(obj_list) <- attrs[["comment"]]
    }
    if (inherits(obj_list, "MizerSim") && isS4(obj_list$params)) {
        obj_list$params <- upgrade_s4_to_s3(obj_list$params)
    }
    obj_list
}

#' S4/S3 compatibility slot access helper
#'
#' @param object An object (S4 or S3).
#' @param name Name of the slot / list element.
#' @param ... Additional arguments.
#' @return The value of the slot or list element.
#' @keywords internal
slot <- function(object, name, ...) {
    if (isS4(object)) {
        object <- upgrade_s4_to_s3(object)
    }
    object[[name]]
}

#' S4/S3 compatibility slot assignment helper
#'
#' @param object An object (S4 or S3).
#' @param name Name of the slot / list element.
#' @param check Logical.
#' @param value Value to assign.
#' @return The modified object.
#' @keywords internal
`slot<-` <- function(object, name, check = TRUE, value) {
    if (isS4(object)) {
        object <- upgrade_s4_to_s3(object)
    }
    object[[name]] <- value
    object
}

#' S4/S3 compatibility .hasSlot helper
#'
#' @param object An object (S4 or S3).
#' @param name Name of the slot / list element.
#' @return Logical.
#' @keywords internal
.hasSlot <- function(object, name) {
    if (isS4(object)) {
        methods::.hasSlot(object, name)
    } else {
        name %in% names(object)
    }
}

#' S4/S3 compatibility slotNames helper
#'
#' @param x An object (S4 or S3) or class name.
#' @return Character vector of slot / element names.
#' @keywords internal
slotNames <- function(x) {
    if (isS4(x)) {
        methods::slotNames(x)
    } else if (is.character(x)) {
        if (identical(x, "MizerParams")) {
            names(emptyParams(c(1, 2), 1))
        } else if (identical(x, "MizerSim")) {
            c("params", "n", "n_pp", "n_other", "effort", "sim_params")
        } else {
            methods::slotNames(x)
        }
    } else if (is.list(x)) {
        names(x)
    } else {
        methods::slotNames(x)
    }
}

#' Validity checker compatibility helper
#'
#' @param object A MizerParams or MizerSim object.
#' @param ... Additional arguments.
#' @return TRUE invisibly if valid, or stops with error.
#' @keywords internal
validObject <- function(object, ...) {
    if (isS4(object)) {
        return(methods::validObject(object, ...))
    }
    if (inherits(object, "MizerParams")) {
        val <- validMizerParams(object)
        if (!isTRUE(val)) {
            stop(paste(val, collapse = "\n"))
        }
        return(invisible(TRUE))
    }
    if (inherits(object, "MizerSim")) {
        val <- valid_MizerSim(object)
        if (!isTRUE(val)) {
            stop(paste(val, collapse = "\n"))
        }
        return(invisible(TRUE))
    }
    invisible(TRUE)
}
