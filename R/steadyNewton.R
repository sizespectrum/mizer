# Direct steady-state solver -------------------------------------------------
#
# `solver = "project"` finds a steady state by running the dynamics until they
# stop changing. That only works for steady states that are dynamically stable:
# at an unstable steady state the time-stepping diverges away from the fixed
# point.
#
# `solver = "newton"` instead solves the discrete steady-state equation
# `F(N) = 0` directly with a Newton-type root finder (via the `nleqslv`
# package). Because it solves the algebraic equation rather than following the
# time evolution, it converges to the steady state irrespective of its dynamic
# stability. See the "Steady-State Solution" section of the
# `vignette("numerical_details")` for the equation that is being solved.
#
# This file holds that solver, `newton_steady_state()`, and the stability
# analysis built on the same machinery. The two exported entry points that call
# the solver, `tuneSteadyState()` and `findSteadyState()`, live in
# `steadyState.R`.

#' Solve the steady-state equation with a Newton-type root finder
#'
#' The engine behind `solver = "newton"` in [tuneSteadyState()] and
#' [findSteadyState()]. It solves the steady-state equation `F(N) = 0` with a
#' globalised Newton iteration from the `nleqslv` package instead of running the
#' dynamics to convergence, so it converges to the steady state even when that
#' steady state is dynamically **unstable**, a case in which the time-stepping
#' solver fails because it diverges away from the fixed point.
#'
#' Two things vary between the two callers, and they are the two arguments that
#' say what is held fixed:
#'
#' * `rdd_const` is the per-species reproduction rate to hold constant, or
#'   `NULL` to let the model's own RDI and RDD functions set it at each iterate.
#'   [tuneSteadyState()] passes the rate of the state it was given;
#'   [findSteadyState()] passes `NULL`.
#' * `resource = "solve"` adds the resource densities to the unknowns and
#'   appends the resource steady-state equation to the residual, so the resource
#'   density and the feeding levels it implies are self-consistent even where
#'   consumers are satiated. The equation is supplied by the
#'   `steady_<resource_dynamics>()` companion and this is what
#'   [findSteadyState()] uses. `resource = "fixed"` leaves the resource out of
#'   the system entirely, holding it at the density in `initial_n_pp`; the
#'   caller is then responsible for rebalancing the capacity afterwards, which
#'   is what [tuneSteadyState()] does.
#'
#' Restoring density-dependent reproduction after a `rdd_const` solve is the
#' caller's job too, so that both solvers share one implementation of it.
#'
#' The consumer densities are solved for in log space, which both keeps them
#' positive and conditions the otherwise badly-scaled system. The unknowns are
#' the size classes in the full potential grid support (running from the egg
#' size up to the grid truncation limit `support_top_idx()`), regardless of
#' whether they carry non-zero density in the supplied initial spectra. A
#' smoothed log-abundance penalty floor is applied to all size classes, which
#' automatically handles zero-density and tail classes smoothly and prevents
#' singular Jacobians. After convergence, size classes that remain at or near
#' the floor are set to zero. This allows the solver to automatically discover
#' the support of the steady state. Newton's method converges from any starting
#' point in the *root's* basin of attraction (which is unrelated to the dynamic
#' stability of the steady state), so a reasonable initial guess should still be
#' supplied in `initialN(params)`.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param resource `"solve"` to make the resource densities unknowns of the
#'   system, `"fixed"` to leave the resource out of it entirely. Its `initial_n` is used as
#'   the starting guess for the iteration.
#' @param effort The fishing effort vector, already validated.
#' @param rdd_const Per-species reproduction rate to hold constant, or `NULL`
#'   to let the reproduction dynamics run.
#' @param resource Either `"solve"` (the resource densities are unknowns of the
#'   system) or `"fixed"` (the resource is held at `initial_n_pp`).
#' @param extinction_floor The relative abundance floor below which a species is
#'   considered extinct. Only used when `rdd_const` is `NULL`.
#' @param verbose Whether to trace the solver iterations to the console.
#' @param solver_tol Tolerance on the residual, passed to
#'   [nleqslv::nleqslv()] as both `ftol` and `xtol`.
#' @param maxit Maximum number of iterations for [nleqslv::nleqslv()].
#' @param jacobian Either `"update"` (the Jacobian is computed once and then
#'   updated cheaply on each iteration) or `"recompute"` (a numerical Jacobian
#'   is computed afresh at every iteration).
#' @param global The globalisation strategy passed to [nleqslv::nleqslv()].
#' @return A list with the new `initial_n` and `initial_n_pp`, the logical
#'   vector of species that went extinct, and the `nleqslv` termination code.
#' @noRd
newton_steady_state <- function(params, effort, rdd_const,
                                resource = c("solve", "fixed"),
                                extinction_floor = 1e-6,
                                verbose = FALSE,
                                solver_tol = 1e-6, maxit = 200,
                                jacobian = c("update", "recompute"),
                                global = "dbldog") {
    resource <- match.arg(resource)
    jacobian <- match.arg(jacobian)
    # nleqslv names these after the Jacobian's origin rather than after what it
    # does; ours says what differs.
    nleqslv_method <- c(update = "Broyden", recompute = "Newton")[[jacobian]]
    if (!requireNamespace("nleqslv", quietly = TRUE)) {
        stop("The Newton solver requires the 'nleqslv' package. ",
             "Install it with install.packages('nleqslv').")
    }

    n_other <- params@initial_n_other
    rates_fns <- projectRateFunctions(params)
    steady_resource <- if (resource == "solve") {
        steady_resource_companion(params)
    } else {
        NULL
    }

    active <- steady_active_set(params, resource = resource)
    if (active$n_fish_active == 0) {
        stop("There is no fish density to solve for: every species is absent ",
             "from `initialN(params)`. The Newton solver holds an absent ",
             "species at zero, so there is nothing left to determine.",
             call. = FALSE)
    }

    # Save initial state for relative floor checks
    N0_initial <- params@initial_n

    # The log-space solve needs a strictly positive, well-scaled start. The
    # second-order schemes can leave isolated zeros inside the support (a
    # negativity-floor artefact), which would make the 1/N-scaled residual
    # overflow. Fill those by log-interpolation from the nonzero neighbours.
    N0 <- positive_initial_guess(params@initial_n, params@w_min_idx,
                                 active$w_top)

    residual_fn <- steady_state_residual(params, rdd_const, n_other, effort,
                                         active,
                                         extinction_floor = extinction_floor,
                                         N0 = N0, rates_fns = rates_fns,
                                         steady_resource = steady_resource)

    x0 <- log(N0[active$mask])
    if (resource == "solve") {
        x0_n_pp <- as.numeric(params@initial_n_pp[active$mask_pp])
        x0_n_pp[x0_n_pp <= 0] <- params@cc_pp[active$mask_pp][x0_n_pp <= 0]
        x0 <- c(x0, log(x0_n_pp))

        initial <- active$unpack(x0)
        resource_eq <- resource_companion_value(params,
                                                n = initial$N,
                                                n_pp = initial$n_pp,
                                                n_other = n_other,
                                                effort = effort,
                                                rates_fns = rates_fns,
                                                steady_resource =
                                                    steady_resource)
        check_resource_equilibrium(resource_eq$n_pp,
                                   active$mask_pp,
                                   steady_resource$name,
                                   where = "at the starting state")
    }

    control <- list(maxit = maxit, ftol = solver_tol, xtol = solver_tol,
                    trace = if (isTRUE(verbose)) 1 else 0)

    sol <- nleqslv::nleqslv(x0, residual_fn, method = nleqslv_method,
                            global = global, control = control)

    unpacked <- active$unpack(sol$x)
    N <- unpacked$N
    n_pp <- unpacked$n_pp
    is_extinct <- rep(FALSE, nrow(N))
    names(is_extinct) <- rownames(N)

    # Check for extinctions if using the relative floor
    if (is.null(rdd_const) && !is.null(extinction_floor) && extinction_floor > 0) {
        extinct_threshold <- extinction_floor * 1.01
        for (i in seq_len(nrow(N))) {
            lo <- params@w_min_idx[i]
            # A species is considered extinct if its recruitment (the egg class)
            # has dropped to the relative abundance floor. Checking any() would
            # produce false positives when the tail naturally shortens.
            if (N0_initial[i, lo] > 0 &&
                N[i, lo] / N0_initial[i, lo] <= extinct_threshold) {
                is_extinct[i] <- TRUE
                N[i, ] <- 0
            }
        }
    }

    # Zero out any tail densities that ended up at or near the structural penalty floor
    for (i in seq_len(nrow(N))) {
        if (!is_extinct[i]) {
            floor_val <- max(N0_initial[i, ]) * 1.1 * 1e-15
            N[i, N[i, ] < floor_val] <- 0
        }
    }

    if (resource == "solve") {
        n_pp <- resource_steady_from_companion(params,
                                               n = N, n_pp = n_pp,
                                               n_other = n_other,
                                               effort = effort,
                                               rates_fns = rates_fns,
                                               steady_resource =
                                                   steady_resource,
                                               mask_pp = active$mask_pp)
    }

    list(n = N, n_pp = n_pp, extinct = is_extinct, absent = active$absent,
         termcd = sol$termcd, message = sol$message)
}

#' The set of size classes solved for by the Newton solver
#'
#' Builds the logical mask of the (species x size) density matrix that the
#' direct solver treats as unknowns. For each species the unknowns run from the
#' egg size `w_min_idx` up to the grid truncation `support_top_idx()`, which is
#' the whole of the size range the dynamics can put density into.
#'
#' The top is the truncation rather than the highest class that currently
#' carries density, because the tail of the solution need not lie under the tail
#' of the starting guess. The growth rate cannot locate the top either: the main
#' reason fish grow beyond `w_repro_max` is diffusion, and the diffusion rate
#' only grows with `w`, never vanishing. Above the truncation the dynamics hold
#' the density at exactly zero, so those classes can never be unknowns; that is
#' what keeps the band of `log(0)` unknowns that would make the Jacobian
#' singular out of the system.
#'
#' Isolated interior zeros (negativity-floor artefacts of the second-order
#' schemes) below the top are kept in the mask and repaired by
#' `positive_initial_guess()`.
#'
#' A species whose row is zero throughout is left out of the mask altogether. It
#' has no density anywhere, so there is nothing to interpolate a starting guess
#' from and `log()` of its row would put `-Inf` into the starting vector. Zero is
#' in any case a steady state of its own equation, so the species is simply held
#' there; recolonisation is not something this solver can discover.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param resource `"solve"` to make the resource densities unknowns of the
#'   system, `"fixed"` to leave them out of it.
#' @return A list with the logical matrix `mask`, the logical vector `absent`
#'   naming the species that were left out, and a function `unpack(x)` that maps
#'   a vector of log-densities (in `mask` order) to the full density matrix.
#' @noRd
steady_active_set <- function(params, resource = c("solve", "fixed")) {
    resource <- match.arg(resource)
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)

    # The grid truncation caps the support: the dynamics hold the density at zero
    # above this, so it can never be an unknown (see support_top_idx()).
    grid_top <- support_top_idx(params)
    n <- params@initial_n
    w_top <- integer(no_sp)
    absent <- rep(FALSE, no_sp)
    names(absent) <- rownames(n)
    mask <- matrix(FALSE, nrow = no_sp, ncol = no_w)
    for (i in seq_len(no_sp)) {
        lo <- params@w_min_idx[i]
        # Always solve up to the grid truncation limit
        w_top[i] <- grid_top[i]
        # A species that is absent to begin with stays absent: it contributes no
        # unknowns, and its row of the solution is left at zero.
        if (all(n[i, ] <= 0)) {
            absent[i] <- TRUE
            next
        }
        mask[i, lo:w_top[i]] <- TRUE
    }
    
    # With the resource held fixed there are no resource unknowns, so the mask
    # is empty and `unpack()` hands back the supplied density untouched.
    if (resource == "solve") {
        mask_pp <- as.logical(params@cc_pp > 0)
    } else {
        mask_pp <- rep(FALSE, length(params@cc_pp))
    }
    n_fish_active <- sum(mask)

    dn <- dimnames(params@initial_n)
    unpack <- function(x) {
        N <- matrix(0, nrow = no_sp, ncol = no_w, dimnames = dn)
        N[mask] <- exp(x[1:n_fish_active])
        n_pp <- params@initial_n_pp
        if (any(mask_pp)) {
            n_pp[mask_pp] <- exp(x[(n_fish_active + 1):length(x)])
        }
        list(N = N, n_pp = n_pp)
    }
    list(mask = mask, w_top = w_top, mask_pp = mask_pp, absent = absent,
         n_fish_active = n_fish_active, unpack = unpack)
}

#' Strictly positive starting guess for the log-space solve
#'
#' The Newton solve uses log-densities as unknowns and scales the residual by
#' `1/N`, so it needs a starting guess that is strictly positive and reasonably
#' scaled on every active size class. The second-order schemes can leave isolated
#' zeros inside the support (a negativity-floor artefact at a reconstructed
#' over/undershoot); a plain `log()` of such a guess is `-Inf`, and flooring it
#' to a tiny constant instead makes the `1/N`-scaled residual overflow. We repair
#' the guess by interpolating `log(N)` linearly across the gaps from the nonzero
#' neighbours (geometric interpolation on the logarithmic size grid), which is
#' both strictly positive and well scaled. Edges with no nonzero neighbour on one
#' side are extrapolated flat.
#'
#' @param N The current density matrix (species x size).
#' @param w_min_idx Per-species egg-size index.
#' @param w_top Per-species support top from `steady_active_set()`.
#' @return A density matrix that is strictly positive on the active range.
#' @noRd
positive_initial_guess <- function(N, w_min_idx, w_top) {
    for (i in seq_len(nrow(N))) {
        rng <- w_min_idx[i]:w_top[i]
        v <- N[i, rng]
        pos <- v > 0
        if (!any(pos)) next
        floor_val <- max(v) * 1e-20
        if (!all(pos)) {
            idx <- seq_along(v)
            v[!pos] <- exp(stats::approx(idx[pos], log(v[pos]),
                                         xout = idx[!pos], rule = 2)$y)
        }
        v[is.na(v) | v <= 0] <- floor_val
        N[i, rng] <- v
    }
    N
}

#' Relative finite-difference step scale
#'
#' Returns the scale that [getStability()] multiplies by `h` to get the
#' finite-difference step for each state variable: the variable's own magnitude
#' where that is nonzero, and otherwise the local scale of the spectrum as
#' supplied by `positive_initial_guess()` (log-interpolated from the nonzero
#' neighbours). Cells whose whole row is zero get no local scale from
#' `positive_initial_guess()`; they fall back to `.Machine$double.eps`, which is
#' harmless because the one-step map is exactly linear about a zero baseline.
#'
#' A step floored at an *absolute* `.Machine$double.eps` instead would be
#' swamped by the rounding error of the order-`N` outputs it perturbs, and the
#' column of the Jacobian would come out as noise (in practice exactly zero).
#'
#' @param x The state variables (a vector of densities).
#' @param local_scale A strictly positive local scale for each entry of `x`, or
#'   zero where none could be determined.
#' @return A strictly positive vector of step scales, the same length as `x`.
#' @noRd
fd_step_scale <- function(x, local_scale) {
    local_scale[!is.finite(local_scale) | local_scale <= 0] <-
        .Machine$double.eps
    pmax(abs(x), local_scale)
}

#' Look up the steady-state companion for resource dynamics
#'
#' Resource dynamics named `foo` opt into the resource-solving Newton branch by
#' supplying a function named `steady_foo()`. It returns the resource abundance
#' that would be steady with the rates held at their supplied values, just as a
#' `balance_foo()` companion supplies parameters that balance those dynamics.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @return A list containing the companion's `name` and resolved function `fun`.
#' @noRd
steady_resource_companion <- function(params) {
    name <- paste0("steady_", params@resource_dynamics)
    fun <- get0(name)
    if (!is.function(fun)) {
        stop("`solver = \"newton\"` cannot solve for the resource with `",
             params@resource_dynamics, "()` because its companion `", name,
             "()` is not available. Define that function or use ",
             "`solver = \"project\"`.", call. = FALSE)
    }
    list(name = name, fun = fun)
}

#' Evaluate a steady-resource companion and the rates it uses
#'
#' Calculates the same complete rate list that the resource dynamics receives
#' during projection, then asks its steady-state companion for the equilibrium
#' implied by those frozen rates. Keeping the rates and the companion result
#' together lets the Newton residual reuse the rates for the consumer equation.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n Consumer densities (species x size).
#' @param n_pp Resource density.
#' @param n_other Abundances of other components.
#' @param effort Fishing effort.
#' @param rates_fns Resolved rate functions from `projectRateFunctions()`.
#' @param steady_resource The result of `steady_resource_companion()`.
#' @return A list containing the equilibrium `n_pp` and the complete `rates`.
#' @noRd
resource_companion_value <- function(params, n, n_pp, n_other, effort,
                                     rates_fns, steady_resource) {
    rates <- rates_fns$Rates(params,
                             n = n, n_pp = n_pp, n_other = n_other, t = 0,
                             effort = effort, rates_fns = rates_fns)
    n_pp_steady <- steady_resource$fun(params,
                                       n = n, n_pp = n_pp,
                                       n_other = n_other, rates = rates, t = 0,
                                       resource_rate = params@rr_pp,
                                       resource_capacity = params@cc_pp)
    if (!is.numeric(n_pp_steady) || !is.null(dim(n_pp_steady)) ||
            length(n_pp_steady) != length(params@w_full)) {
        stop("`", steady_resource$name, "()` must return a numeric vector of ",
             "length ", length(params@w_full), ", one value for each resource ",
             "size class.", call. = FALSE)
    }
    list(n_pp = as.numeric(n_pp_steady), rates = rates)
}

#' Require the positive resource branch used by the Newton solver
#'
#' Resource densities are log-space unknowns and are therefore strictly
#' positive. A companion result at or below zero belongs to a boundary branch
#' that this solver does not yet represent.
#'
#' @param n_pp_steady The equilibrium returned by a resource companion.
#' @param mask_pp Resource size classes included in the Newton system.
#' @param companion_name Name of the resource companion.
#' @param where Description of the state at which it was evaluated.
#' @return `NULL` invisibly, or an error.
#' @noRd
check_resource_equilibrium <- function(n_pp_steady, mask_pp,
                                       companion_name, where) {
    bad <- mask_pp & (!is.finite(n_pp_steady) | n_pp_steady <= 0)
    if (!any(bad)) {
        return(invisible(NULL))
    }
    stop("`", companion_name, "()` does not give a finite, positive resource ",
         "equilibrium in ", sum(bad), " size ",
         if (sum(bad) == 1) "class " else "classes ", where, ". The Newton ",
         "solver currently supports only positive resource equilibria; use ",
         "`solver = \"project\"` when resource size classes can be depleted ",
         "to zero.", call. = FALSE)
}

#' Re-equilibrate the resource through its steady-state companion
#'
#' Repeats the frozen-rate equilibrium calculation because consumer satiation
#' makes resource mortality depend indirectly on resource abundance. This is
#' the dynamics-agnostic replacement for the former semichemostat equation.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n Consumer densities (species x size).
#' @param n_pp Resource density from the coupled Newton solve.
#' @param n_other Abundances of other components.
#' @param effort Fishing effort.
#' @param rates_fns Resolved rate functions from `projectRateFunctions()`.
#' @param steady_resource The result of `steady_resource_companion()`.
#' @param mask_pp Resource size classes included in the Newton system.
#' @return The steady-state resource number density vector.
#' @noRd
resource_steady_from_companion <- function(params, n, n_pp, n_other, effort,
                                           rates_fns, steady_resource,
                                           mask_pp) {
    for (i in 1:8) {
        resource_eq <- resource_companion_value(params,
                                                n = n, n_pp = n_pp,
                                                n_other = n_other,
                                                effort = effort,
                                                rates_fns = rates_fns,
                                                steady_resource =
                                                    steady_resource)
        check_resource_equilibrium(resource_eq$n_pp,
                                   mask_pp, steady_resource$name,
                                   where = "after the Newton solve")
        n_pp <- resource_eq$n_pp
    }
    n_pp
}

#' Unscaled residual of the discrete steady-state equation
#'
#' Evaluates the growth, mortality and diffusion rates at the given state,
#' assembles the steady-state tridiagonal coefficients with
#' `get_transport_coefs()` at `dt = 1` and returns
#' \deqn{F_j = a_j N_{j-1} + b_j N_j + c_j N_{j+1} - S_j.}
#' Since `S = N + recruitment source`, the `+N` cancels the backward-Euler
#' `N^t` term, leaving exactly the steady-state equation
#' \eqn{\tilde A N_{j-1} + \tilde B N_j + \tilde C N_{j+1} - \tilde S = 0}.
#'
#' Because the backward-Euler coefficients are `A = I - dt L`, this residual is
#' \eqn{-dt\,dN/dt} exactly, so at `dt = 1` its negative is the instantaneous
#' rate of change of the density. That is what makes it usable both as the
#' function whose root the Newton solver seeks and as the diagnostic
#' [getSteadyResidual()] reports.
#'
#' This is the single implementation of the transport residual; both callers go
#' through it.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n Consumer densities (species x size) to evaluate the residual at.
#' @param n_pp Resource density to evaluate the rates at.
#' @param n_other Abundances of other components.
#' @param effort The fishing effort vector.
#' @param rdd Per-species reproduction rate to use. If `NULL`, it is calculated
#'   from the model's own RDI and RDD functions at the given state.
#' @param rates_fns The rate functions, from `projectRateFunctions()`. Passed in
#'   so that a caller evaluating the residual many times need only look them up
#'   once.
#' @param flux_limiter The flux limiter scheme, from `flux_limiter_scheme()`.
#' @param rates The rates at this state, as returned by `mizer_rates_subset()`
#'   with targets `EGrowth`, `Mort` and `Diffusion`. Pass them in when the
#'   caller has computed them already; otherwise they are computed here.
#' @return The unscaled residual matrix (species x size).
#' @noRd
consumer_residual <- function(params, n, n_pp, n_other, effort, rdd = NULL,
                              rates_fns = projectRateFunctions(params),
                              flux_limiter = flux_limiter_scheme(params),
                              rates = NULL) {
    no_w <- length(params@w)

    r <- rates
    if (is.null(r)) {
        r <- mizer_rates_subset(params, n = n, n_pp = n_pp, n_other = n_other,
                                t = 0, effort = effort, rates_fns = rates_fns,
                                targets = c("EGrowth", "Mort", "Diffusion"))
    }

    if (is.null(rdd)) {
        rdd <- state_rdd(params, n = n, n_pp = n_pp, n_other = n_other,
                         rates = r)
    }

    coefs <- get_transport_coefs(params, n = n, g = r$e_growth,
                                 mu = r$mort, dt = 1,
                                 recruitment_flux = rdd,
                                 d = r$diffusion,
                                 flux_limiter = flux_limiter)

    Nm <- cbind(0, n[, -no_w, drop = FALSE])   # N_{j-1}
    Np <- cbind(n[, -1, drop = FALSE], 0)      # N_{j+1}
    coefs$a * Nm + coefs$b * n + coefs$c * Np - coefs$S
}

#' Reproduction rate at a given state
#'
#' Runs the model's own RDI and RDD functions at the state `(n, n_pp, n_other)`,
#' going through the extension dispatch when the model carries extensions. This
#' is the reproduction rate that the dynamics would use at that state, as
#' opposed to a rate held fixed at its steady-state value. It is shared by
#' `consumer_residual()` and the stability analyses, which have to agree on it.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param n Consumer densities (species x size).
#' @param n_pp Resource density.
#' @param n_other Abundances of other components.
#' @param rates The rates at this state, at least `e_repro`, `e_growth`, `mort`
#'   and `diffusion`, as returned by `mizer_rates_subset()`.
#' @return The per-species density-dependent reproduction rate.
#' @noRd
state_rdd <- function(params, n, n_pp, n_other, rates) {
    if (usesExtensionDispatch(params)) {
        rdi <- projectRDI(params, n = n, n_pp = n_pp, n_other = n_other, t = 0,
                          e_repro = rates$e_repro, e_growth = rates$e_growth,
                          mort = rates$mort, diffusion = rates$diffusion)
        projectRDD(params, rdi = rdi, species_params = params@species_params,
                   t = 0)
    } else {
        f_rdi <- get(params@rates_funcs$RDI)
        rdi <- f_rdi(params, n = n, n_pp = n_pp, n_other = n_other, t = 0,
                     e_repro = rates$e_repro, e_growth = rates$e_growth,
                     mort = rates$mort, diffusion = rates$diffusion)
        f_rdd <- get(params@rates_funcs$RDD)
        f_rdd(rdi = rdi, species_params = params@species_params,
              params = params, t = 0)
    }
}

#' Residual of the discrete steady-state equation
#'
#' Returns a closure `f(x)` suitable for [nleqslv::nleqslv()]. The argument `x`
#' is the vector of log-densities of the active size classes (see
#' `steady_active_set()`). The closure rebuilds the full density matrix,
#' calls `consumer_residual()` for the consumer steady-state residual
#' \eqn{F_j = a_j N_{j-1} + b_j N_j + c_j N_{j+1} - S_j}, to which it adds the
#' scaling, masking and penalty floor that the root finder needs. When the
#' resource is part of the solve, the closure obtains its frozen-rate
#' equilibrium from the `steady_<resource_dynamics>()` companion and appends
#' the relative difference between that equilibrium and the current resource.
#' The residual is divided by `N`, turning it into a per-capita rate of change
#' that is dimensionless and O(1) across the many orders of magnitude spanned by
#' the densities — the natural scaling to pair with the log-space unknowns.
#'
#' The flux-limiter weight is recomputed fresh on every evaluation, so the van
#' Leer / centred nonlinearity is part of the root and the solution matches what
#' [project()] converges to on the same `params`.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param rdd_const Per-species reproduction rate held constant during the solve.
#' @param n_other Abundances of other components (held constant).
#' @param effort The fishing effort vector.
#' @param active The active-set list from `steady_active_set()`.
#' @param rates_fns Resolved rate functions from `projectRateFunctions()`.
#' @param steady_resource The result of `steady_resource_companion()`, or `NULL`
#'   when the resource is fixed.
#' @return A function of the packed log-density vector returning the packed
#'   scaled residual.
#' @noRd
steady_state_residual <- function(params, rdd_const, n_other, effort, active,
                                  extinction_floor = 1e-6, N0 = NULL,
                                  rates_fns = projectRateFunctions(params),
                                  steady_resource = NULL) {
    no_w <- length(params@w)
    mask <- active$mask
    mask_pp <- active$mask_pp
    n_fish_active <- active$n_fish_active
    flux_limiter <- flux_limiter_scheme(params)

    if (is.null(N0)) {
        x0_initial <- params@initial_n
    } else {
        x0_initial <- N0
    }

    # A floor is always active to handle zero-abundance tail classes smoothly
    # and prevent singular Jacobians.
    x0 <- log(x0_initial[mask])
    support_floor <- matrix(0, nrow = nrow(x0_initial), ncol = no_w)
    for (i in seq_len(nrow(x0_initial))) {
        support_floor[i, ] <- log(max(x0_initial[i, ])) + log(1e-15)
    }
    x_floor <- support_floor[mask]

    if (is.null(rdd_const) && !is.null(extinction_floor) &&
            extinction_floor > 0) {
        ext_floor <- x0 + log(extinction_floor)
        x_floor <- pmax(x_floor, ext_floor)
    }

    function(x) {
        unpacked <- active$unpack(x)
        N <- unpacked$N
        n_pp <- unpacked$n_pp

        resource_eq <- NULL
        rates <- NULL
        if (any(mask_pp)) {
            resource_eq <- resource_companion_value(params,
                                                    n = N, n_pp = n_pp,
                                                    n_other = n_other,
                                                    effort = effort,
                                                    rates_fns = rates_fns,
                                                    steady_resource =
                                                        steady_resource)
            if (any(!is.finite(resource_eq$n_pp))) {
                stop("`", steady_resource$name, "()` returned non-finite ",
                     "resource equilibria during the Newton iteration.",
                     call. = FALSE)
            }
            rates <- resource_eq$rates
        }

        res <- consumer_residual(params, n = N, n_pp = n_pp, n_other = n_other,
                                 effort = effort, rdd = rdd_const,
                                 rates_fns = rates_fns,
                                 flux_limiter = flux_limiter, rates = rates)

        # Scale to a per-capita rate of change (N > 0 on the active set).
        if (is.null(rdd_const)) {
            r_ss <- (res / x0_initial)[mask]
        } else {
            r_ss <- (res / N)[mask]
        }

        if (!is.null(x_floor)) {
            y <- x[1:n_fish_active] - x_floor
            delta <- 0.001
            penalty <- 0.5 * (y - sqrt(y^2 + delta^2))
            r_ss <- r_ss + penalty
        }

        # Resource residual, only when the resource is part of the system.
        if (!any(mask_pp)) {
            return(r_ss)
        }
        r_pp_ss <- ((resource_eq$n_pp - n_pp) / n_pp)[mask_pp]

        c(r_ss, as.numeric(r_pp_ss))
    }
}

# Stability analysis ----------------------------------------------------------

#' Shared scaffolding for the stability analyses
#'
#' [getStability()] and [getDiscreteStability()] differ only in the map they
#' linearise: the continuous rates of change in the one case, mizer's numerical
#' one-step map in the other. Everything around that map is the same, and this
#' function assembles it: the validated `params`, the active set, the packed
#' state vector with its finite-difference step scales, and the closures that
#' unpack a packed state, name a cell and check a result for non-finite values.
#'
#' @param params A \linkS4class{MizerParams} object.
#' @param effort The fishing effort to use.
#' @param fn_name The name of the calling function, for error messages.
#' @param map_name What the linearised map is called, for error messages.
#' @return A list holding the validated `params`, the state and the closures
#'   described above.
#' @noRd
stability_context <- function(params,
                              effort = params@initial_effort,
                              fn_name = "getStability",
                              map_name = "the rates of change") {
    params <- validParams(params)
    effort <- validEffortVector(effort, params = params)
    params@initial_effort <- effort

    # The whole calculation linearises the dynamics *at* the stored state. If
    # that state is not a fixed point, the eigenvalues describe the
    # neighbourhood of a point the model is not sitting at, and the verdict on
    # stability means nothing.
    warn_if_not_steady(params, paste(
        "The eigenvalues below describe the dynamics near a fixed point, so on",
        "a state that is not one they are not meaningful. Use",
        "`findSteadyState()`",
        "first."))

    # The Jacobian below has a row for every fish cell and every resource cell,
    # and none for any other component: they enter as fixed inputs. A component
    # that responds to the fish on a comparable timescale is then missing from
    # the spectrum, which can change the verdict.
    warn_other_components_fixed(params, paste(
        "The eigenvalues therefore describe the consumer-resource subsystem",
        "with", if (length(params@other_dynamics) > 1) "them" else "it",
        "frozen, which is the whole story only if the component responds to",
        "the fish much faster or much more slowly than the fish respond to",
        "each other."))

    n_other      <- params@initial_n_other
    active       <- steady_active_set(params)
    active_idx   <- which(active$mask)
    flux_limiter <- flux_limiter_scheme(params)
    rates_fns    <- projectRateFunctions(params)

    N_ss   <- params@initial_n
    npp_ss <- params@initial_n_pp
    N_vec  <- N_ss[active$mask]
    n_fish_active <- length(N_vec)

    # Finite-difference step sizes. The step is relative, `h * N`, so a cell at
    # exactly zero would get a zero step: the difference quotient then drops
    # the cell from the Jacobian and puts a spurious zero eigenvalue in place of
    # its true decay rate. Such cells do occur: an isolated negativity-floor
    # artefact of the second-order schemes, or a tail class that the solver
    # zeroed at its structural floor. We therefore floor the step at the local
    # scale of the spectrum, log-interpolated across the gaps from the nonzero
    # neighbours exactly as the Newton solve does. A row that is zero throughout
    # (an extinct species) has no local scale either, but there the dynamics are
    # exactly linear about a zero baseline, so any step recovers the derivative
    # and the absolute floor is kept.
    N_scale <- fd_step_scale(N_vec,
                             positive_initial_guess(N_ss, params@w_min_idx,
                                                    active$w_top)[active$mask])

    npp_vec <- as.numeric(npp_ss)
    n_npp   <- length(npp_vec)
    # Step scale for the resource, floored the same way as for the fish.
    # The resource carries structural zeros above `w_pp_cutoff`, where the
    # capacity vanishes; log-interpolation with flat extrapolation gives
    # them the scale of the last nonzero class.
    npp_local <- positive_initial_guess(matrix(npp_vec, nrow = 1),
                                        w_min_idx = 1L, w_top = n_npp)
    npp_scale <- fd_step_scale(npp_vec, as.numeric(npp_local))
    x0    <- c(N_vec, npp_vec)
    scale <- c(N_scale, npp_scale)

    # Packed state -> the arrays the rate functions want.
    unpack <- function(x) {
        N <- N_ss
        N[active_idx] <- x[seq_len(n_fish_active)]
        n_pp <- npp_ss
        n_pp[] <- x[n_fish_active + seq_len(n_npp)]
        list(N = N, n_pp = n_pp)
    }

    # A custom rate function can be non-finite at a state that the differencing
    # visits. Fail loudly, naming the cell, rather than letting a NaN travel on
    # into eigen() and come back as a meaningless spectrum.
    check <- function(x, where) {
        if (!all(is.finite(x))) {
            stop(fn_name, "(): ", map_name, " returned non-finite values ",
                 where, ". Every state it evaluates satisfies N >= 0, so this ",
                 "points to a rate function that is not finite there.",
                 call. = FALSE)
        }
    }

    describe <- function(k) {
        if (k <= n_fish_active) {
            rc <- arrayInd(active_idx[k], dim(N_ss))
            paste0("when perturbing species ",
                   params@species_params$species[rc[1]], " at w = ",
                   signif(params@w[rc[2]], 3))
        } else {
            paste0("when perturbing the resource at w = ",
                   signif(params@w_full[k - n_fish_active], 3))
        }
    }

    list(params = params, effort = effort,
         n_other = n_other, active = active, active_idx = active_idx,
         n_fish_active = n_fish_active, n_npp = n_npp,
         rates_fns = rates_fns, flux_limiter = flux_limiter,
         N_ss = N_ss, npp_ss = npp_ss, x0 = x0, scale = scale,
         unpack = unpack, check = check, describe = describe)
}

#' Finite-difference Jacobian of a map on the packed state
#'
#' Differences `f` column by column at the steady state, with a multiplicative
#' step `h * scale[k]`.
#'
#' Every state at which `f` is evaluated satisfies \eqn{N \ge 0}: where a
#' centred step would push a cell negative — which can only happen for a cell at
#' (or below) the floor described in `stability_context()` — the column is
#' differenced forwards from the unperturbed state instead. Such columns are
#' first order in `h` rather than second.
#'
#' @param ctx The context from `stability_context()`.
#' @param f The map, a function of the packed state vector returning a vector of
#'   the same length.
#' @param h Relative step size.
#' @return The Jacobian matrix of `f` at the steady state.
#' @noRd
stability_jacobian <- function(ctx, f, h) {
    x <- ctx$x0
    n <- length(x)
    L <- matrix(0, nrow = n, ncol = n)

    # Baseline for the one-sided columns, shared by all of them. It has to be
    # evaluated rather than assumed: the state is a fixed point only to the
    # tolerance of whatever produced it, and that residual would otherwise be
    # divided by eps into every one-sided column.
    base <- f(x)
    ctx$check(base, "at the unperturbed steady state")

    for (k in seq_len(n)) {
        eps <- h * ctx$scale[k]
        x_plus <- x
        x_plus[k] <- x[k] + eps
        if (x[k] - eps < 0) {
            # A centred step would leave the physical cone N >= 0. Difference
            # forwards instead: at the boundary of the cone the one-sided
            # derivative is the right object anyway, because the dynamics never
            # visit the states a centred step would sample. This keeps every
            # rate function evaluation inside N >= 0, so extensions need not be
            # defined at negative abundances.
            L[, k] <- (f(x_plus) - base) / eps
        } else {
            x_minus <- x
            x_minus[k] <- x[k] - eps
            L[, k] <- (f(x_plus) - f(x_minus)) / (2 * eps)
        }
        ctx$check(L[, k], ctx$describe(k))
    }
    L
}

#' The continuous rates of change as a map on the packed state
#'
#' Returns the function that [getStability()] differentiates: the instantaneous
#' rate of change of every state variable. For the consumers this is
#' `-consumer_residual()`, which assembles exactly the semi-discretised
#' \eqn{dN/dt} with the spatial scheme configured in `params`. For the resource
#' it is the analytic semichemostat derivative; any other resource dynamics is
#' only available as a one-step map, so there it is differenced over a short
#' step, as `steady_rates()` does.
#'
#' @param ctx The context from `stability_context()`.
#' @param dt_resource Step length for differencing a resource dynamics function
#'   that is not the semichemostat.
#' @return A function of the packed state vector.
#' @noRd
stability_rhs <- function(ctx, dt_resource = 1e-4) {
    params  <- ctx$params
    targets <- c("EGrowth", "Mort", "Diffusion", "ResourceMort")
    semichemostat <- params@resource_dynamics == "resource_semichemostat"
    resource_dyn  <- if (semichemostat) NULL else get(params@resource_dynamics)

    function(x) {
        st <- ctx$unpack(x)
        r <- mizer_rates_subset(params, n = st$N, n_pp = st$n_pp,
                                n_other = ctx$n_other, t = 0,
                                effort = ctx$effort,
                                rates_fns = ctx$rates_fns, targets = targets)
        # `rdd = NULL` makes `consumer_residual()` evaluate the model's own
        # reproduction function at the perturbed state, so the reproduction
        # feedback is part of the Jacobian just as it is part of `project()`.
        dNdt <- -consumer_residual(params, n = st$N, n_pp = st$n_pp,
                                   n_other = ctx$n_other, effort = ctx$effort,
                                   rdd = NULL,
                                   rates_fns = ctx$rates_fns,
                                   flux_limiter = ctx$flux_limiter,
                                   rates = r)
        if (semichemostat) {
            # dn_pp/dt = rr (cc - n_pp) - mu_R n_pp. Where both the
            # replenishment rate and the mortality vanish this is zero, so
            # such a class simply stays put, as resource_semichemostat()
            # also arranges.
            dnpp <- params@rr_pp * params@cc_pp -
                (params@rr_pp + r$resource_mort) * st$n_pp
        } else {
            npp_new <- resource_dyn(params, n = st$N, n_pp = st$n_pp,
                                    n_other = ctx$n_other, rates = r,
                                    t = 0, dt = dt_resource,
                                    resource_rate = params@rr_pp,
                                    resource_capacity = params@cc_pp)
            dnpp <- (npp_new - st$n_pp) / dt_resource
        }
        c(dNdt[ctx$active$mask], as.numeric(dnpp))
    }
}

#' Mizer's numerical one-step map on the packed state
#'
#' Returns the function that [getDiscreteStability()] differentiates: one step
#' of length `dt` of the dynamics as [project()] takes it with
#' `method = "euler"`, using the same `project_n_loop()` C++ Thomas solver, the
#' same transport coefficients and the model's own `resource_dynamics`.
#'
#' It has to be that map and not one that merely agrees with it to first order
#' in `dt`, because the whole purpose of the function is to say what mizer's
#' solver does at a particular `dt`. An earlier version stepped the resource
#' with backward Euler so that the fast resource modes could not land on a
#' discrete eigenvalue of exactly zero. There is nothing wrong with such an
#' eigenvalue: it is the correct statement that the mode relaxes inside a single
#' step, and it can never be the largest in modulus, so it never touches the
#' spectral radius.
#'
#' @param ctx The context from `stability_context()`.
#' @param dt The step length.
#' @return A function of the packed state vector.
#' @noRd
stability_step <- function(ctx, dt) {
    params  <- ctx$params
    w_top   <- support_top_idx(params)
    targets <- c("EGrowth", "Mort", "Diffusion", "ResourceMort")
    resource_dyn <- get(params@resource_dynamics)

    function(x) {
        st <- ctx$unpack(x)
        r <- mizer_rates_subset(params, n = st$N, n_pp = st$n_pp,
                                n_other = ctx$n_other, t = 0,
                                effort = ctx$effort,
                                rates_fns = ctx$rates_fns, targets = targets)
        rdd <- state_rdd(params, n = st$N, n_pp = st$n_pp,
                         n_other = ctx$n_other, rates = r)
        coefs <- get_transport_coefs(params, n = st$N, g = r$e_growth,
                                     mu = r$mort, dt = dt,
                                     recruitment_flux = rdd,
                                     d = r$diffusion,
                                     flux_limiter = ctx$flux_limiter)
        N_out <- project_n_loop(st$N, coefs$a, coefs$b, coefs$c, coefs$S,
                                params@w_min_idx)
        npp_out <- resource_dyn(params, n = st$N, n_pp = st$n_pp,
                                n_other = ctx$n_other, rates = r,
                                t = 0, dt = dt,
                                resource_rate = params@rr_pp,
                                resource_capacity = params@cc_pp)
        c(zero_above_support(N_out, w_top)[ctx$active$mask],
          as.numeric(npp_out))
    }
}

#' Reshape the leading eigenvectors back into abundance space
#'
#' @param evecs The eigenvector matrix, columns already sorted.
#' @param ctx The context from `stability_context()`.
#' @param cols The columns of `evecs` to reshape. Defaults to the first two.
#' @return A list with `$fish`, a complex `(n_species, n_sizes, length(cols))`
#'   array, and `$resource`, a complex `(n_w_full, length(cols))` matrix.
#' @noRd
stability_eigenvectors <- function(evecs, ctx,
                                   cols = seq_len(min(2L, ncol(evecs)))) {
    N_ss  <- ctx$N_ss
    no_sp <- nrow(N_ss)
    no_w  <- ncol(N_ss)
    n_out <- length(cols)

    fish <- array(0 + 0i, dim = c(no_sp, no_w, n_out),
                  dimnames = c(dimnames(N_ss), list(NULL)))
    resource <- matrix(0 + 0i, nrow = ctx$n_npp, ncol = n_out)

    # Fish abundances and resource densities differ by many orders of
    # magnitude, so a normalisation by absolute modulus would be set entirely
    # by the resource block and leave the fish block at ~1e-13. Scale instead
    # by the largest perturbation *relative to the local steady state*, which
    # is the quantity the two blocks have in common. `ctx$x0` is the packed
    # steady state, in the same row order as the eigenvectors.
    x_ss <- ctx$x0
    pos  <- x_ss > 0

    for (k in seq_len(n_out)) {
        # One scalar for the whole state vector, so that the relative
        # amplitude and phase between the fish block and the resource block
        # are those of the mode. Normalising the two blocks separately would
        # rescale them by different factors and destroy exactly the
        # information that makes the resource component worth returning.
        v <- evecs[, cols[k]]
        max_rel <- if (any(pos)) max(Mod(v[pos]) / x_ss[pos]) else 0
        if (max_rel > 0) v <- v / max_rel

        M <- matrix(0 + 0i, nrow = no_sp, ncol = no_w)
        M[ctx$active_idx] <- v[seq_len(ctx$n_fish_active)]
        fish[, , k] <- M
        resource[, k] <- v[ctx$n_fish_active + seq_len(ctx$n_npp)]
    }
    rownames(resource) <- names(ctx$npp_ss)
    list(fish = fish, resource = resource)
}

#' Drop the trailing mode dimension of a `stability_eigenvectors()` result
#'
#' @param ev A `stability_eigenvectors()` list holding a single mode.
#' @return A list with `$fish`, a complex `(n_species, n_sizes)` matrix, and
#'   `$resource`, a complex vector of length `n_w_full`.
#' @noRd
single_eigenvector <- function(ev) {
    list(fish = ev$fish[, , 1], resource = ev$resource[, 1])
}

#' Analyse the dynamic stability of a mizer steady state
#'
#' `r lifecycle::badge("experimental")`
#' Computes the eigenvalues of the linearised dynamics at the steady state
#' stored in `params@initial_n`. These eigenvalues determine whether the steady
#' state is dynamically stable and, where the spectrum contains a complex pair,
#' the period at which the model oscillates.
#'
#' ## Mathematical background
#'
#' Mizer discretises the size axis but not time: on the size grid the model is a
#' system of ordinary differential equations
#' \deqn{\frac{dN}{dt} = F(N, n_{pp}),}
#' where \eqn{F} collects the divergence of the growth flux, the mortality sink
#' and the reproductive influx at the egg size, assembled with the spatial
#' scheme configured via [second_order_w()]. `getStability()` differentiates
#' \eqn{F} directly, by centred finite differences in each state variable, and
#' returns the eigenvalues \eqn{\lambda_i} of the resulting Jacobian
#' \eqn{J = \partial F/\partial N}. The steady state is **stable** when all of
#' them satisfy \eqn{\text{Re}(\lambda_i) < 0} and **unstable** when at least
#' one exceeds 0.
#'
#' No time step enters this calculation. The eigenvalues are a property of the
#' model, not of any solver: they describe the continuous-time dynamics of the
#' semi-discretised model, and are what a simulation with a small enough time
#' step converges to. The stability of the numerical step itself is a separate
#' question, answered by [getDiscreteStability()].
#'
#' A complex-conjugate pair \eqn{\lambda = \sigma \pm i\omega} is an
#' oscillatory mode: a perturbation along it rings with period
#' \deqn{T = \frac{2\pi}{|\omega|} \text{ years,}}
#' growing or decaying as \eqn{e^{\sigma t}}. The pair with the largest
#' \eqn{\sigma} is returned as `leading_oscillatory_eigenvalue`, with its
#' period and eigenvector.
#'
#' That is a statement about the mode, not about a bifurcation. A **Hopf
#' bifurcation** is the event of such a pair *crossing* the imaginary axis, and
#' a single spectrum cannot show a crossing: the leading oscillatory mode of a
#' comfortably stable model can sit far to the left, ringing only as a transient
#' on the way back to the fixed point. Establishing a Hopf bifurcation means
#' watching \eqn{\sigma} pass through zero as a parameter is varied, which is
#' what [scanModel()] is for. Only then is \eqn{T} the period of an emerging
#' limit cycle; otherwise it is the period of a damped (or growing) oscillation.
#'
#' ## What is in the Jacobian
#'
#' The resource is a state variable of the system like any other: fish and
#' resource cells are perturbed independently, giving the full coupled Jacobian.
#' Its eigenvalues include both the slow fish modes and a cluster of fast
#' resource-relaxation modes, at \eqn{\lambda \approx -(r_{pp} + \mu_R)}. Any
#' resource dynamics function is supported: the semichemostat derivative is
#' written down analytically, and anything else is differenced over a short step.
#'
#' Components registered with [setComponent()] are *not* state variables here.
#' They are held at their stored values while the fish and the resource are
#' perturbed, so the spectrum is that of the consumer-resource subsystem with
#' the components frozen. This is exact when a component is a fixed input, and a
#' good approximation when it is much faster or much slower than the fish, but
#' it is not the full model, and mizer says so with a warning when it meets one.
#' Giving extension components an explicit residual and Jacobian is the work
#' that would lift this restriction.
#'
#' Reproduction is a state-dependent rate like any other: the reproduction
#' function stored in `params@rates_funcs$RDD` is evaluated at each perturbed
#' state, so the feedback from the spectra back onto the influx of eggs is part
#' of the Jacobian, exactly as it is part of [project()]. There is no option to
#' pin the reproduction rate at its value at the fixed point. A model in which
#' reproduction really is constant expresses that as a model: with
#' `rates_funcs$RDD = "constantRDD"` the derivative of the reproduction rate is
#' zero and the pinned Jacobian is what the analysis returns.
#'
#' This is why the stability of a steady state depends on the reproduction
#' parameters even though the steady state itself does not.
#' [setBevertonHolt()] moves along a family of `erepro`/`R_max` pairs that all
#' leave the same fixed point, but they do not all leave the same dynamics: at a
#' [reproduction_level()] near 1 the reproduction rate barely responds to the
#' energy invested in it, approaching the constant-reproduction case, while at a
#' level near 0 it follows that energy proportionally. The two ends can differ
#' in their verdict, so the analysis has to read the model rather than take an
#' argument.
#'
#' ## Numerical details
#'
#' The Jacobian is computed numerically using a multiplicative (relative)
#' finite-difference step \eqn{h \cdot N^*}. Where a cell sits at exactly zero
#' and so has no scale of its own, the step is floored at the local scale of the
#' spectrum, interpolated from the nonzero neighbours, so that the cell still
#' gets a resolved derivative rather than a column of rounding error.
#'
#' Every state at which the rates are evaluated satisfies \eqn{N \ge 0}: where
#' a centred step would push a cell negative — which can only happen for a cell
#' at (or below) the floor described above — the column is differenced forwards
#' from the unperturbed state instead. At the boundary of the physical cone the
#' one-sided derivative is the appropriate object anyway, since the dynamics
#' never visit the states a centred step would sample. A rate function
#' registered with [setRateFunction()] therefore never has to be defined at
#' negative abundances. Such columns are first order in `h` rather than second,
#' so they respond slightly more to a change of `h` than the rest.
#'
#' @param params A \linkS4class{MizerParams} object whose `initial_n` holds the
#'   steady state to analyse. Typically the output of [findSteadyState()].
#' @param effort The fishing effort to use. By default the initial effort
#'   stored in `params`.
#' @param h Relative step size for centred finite differences. Default `1e-4`.
#'   The result should not depend on this choice. If it does, the dynamics are
#'   not smooth at the state being analysed — see the section below.
#' @return A named list with the following components:
#'   \describe{
#'     \item{`eigenvalues`}{Complex vector of the continuous-time eigenvalues
#'       \eqn{\lambda_i}, sorted by decreasing real part.}
#'     \item{`max_real_part`}{The largest real part of the eigenvalues:
#'       \eqn{\max_i \text{Re}(\lambda_i)}. Greater than 0 means unstable.}
#'     \item{`stable`}{Logical: `TRUE` when `max_real_part < 0`.}
#'     \item{`dominant_period`}{The period (in years) of the dominant
#'       eigenvalue: `2*pi / abs(Im(lambda_1))`. `Inf` for a real
#'       dominant eigenvalue (monotone dynamics).}
#'     \item{`oscillation_period`}{Period (in years) of the oscillatory mode
#'       below, \eqn{2\pi/|\omega|}; `NULL` when no complex eigenvalue exists.
#'       It is the period at which the model rings, and only at a Hopf
#'       bifurcation — where the real part is zero — the period of a limit
#'       cycle.}
#'     \item{`leading_oscillatory_eigenvalue`}{The complex eigenvalue with the
#'       largest real part, or `NULL` when there is none. Its real part is the
#'       rate at which that oscillation grows, so a strongly negative one means
#'       the ringing is a transient, not a cycle the model settles onto.}
#'     \item{`leading_oscillatory_eigenvector`}{Its eigenvector, as a list with `$fish`, a
#'       complex `(n_species, n_sizes)` matrix, and `$resource`, a complex
#'       vector of length `n_w_full`. This is the mode [getOscillationModeSim()]
#'       draws, and it is not in general one of `leading_eigenvectors`: the
#'       dominant mode of the system can be real while the dominant
#'       *oscillatory* mode is well down the spectrum.}
#'     \item{`n_active`}{Dimension of the Jacobian: the number of active fish
#'       cells plus all resource cells.}
#'     \item{`leading_eigenvectors`}{The eigenvectors of the two eigenvalues
#'       with the largest real part, reshaped back into the state space: a list
#'       with `$fish`, a complex array of shape `(n_species, n_sizes, 2)` with
#'       the same species and size dimnames as `params@initial_n`, and
#'       `$resource`, a complex matrix of shape `(n_w_full, 2)`.
#'       Each eigenvector is normalised by a single scalar covering both
#'       blocks, so that the relative amplitude and phase between fish and
#'       resource are those of the mode. The scalar is chosen so that the
#'       largest perturbation *relative to the steady state*,
#'       \eqn{|v_i| / x^*_i}, is 1 somewhere in the state: an absolute
#'       normalisation would be set entirely by the resource, whose densities
#'       dwarf the fish abundances. `Mod(fish[, , 1]) / initialN(params)` is
#'       therefore the relative amplitude pattern, peaking at 1 in whichever
#'       cell swings hardest. The real and imaginary parts of eigenvector 1
#'       span the two-dimensional oscillation plane of the dominant mode.}
#'     \item{`params`}{The validated `params` object the analysis was made at.}
#'   }
#' @section Requires smooth dynamics:
#' The finite-difference Jacobian is only meaningful if the rates of change are
#' differentiable at \eqn{N^*}. A custom rate function registered with
#' [setRateFunction()] that jumps as a function of the abundances breaks this in
#' two ways. If the state sits on the switching threshold, some perturbations
#' straddle it and pick up the jump, and the reported eigenvalues then vary
#' wildly with `h`. If the state is near but not on the threshold, no
#' perturbation crosses it, and the function silently returns the stability of
#' the single branch the state happens to lie on — which can read as `stable`
#' for a model whose simulations never settle.
#'
#' Re-running with a different `h` is the cheapest check: if the answer moves,
#' do not trust it. See [Discontinuous rate
#' functions](https://sizespectrum.org/mizer/articles/discontinuous_rates.html).
#'
#' @seealso [findSteadyState()], [getDiscreteStability()], [getOscillationModeSim()]
#' @export
getStability <- function(params,
                         effort = params@initial_effort,
                         h = 1e-4) {
    assert_that(is.number(h), is.finite(h), h > 0)
    ctx <- stability_context(params,
                             effort = effort,
                             fn_name = "getStability",
                             map_name = "the rates of change")
    J <- stability_jacobian(ctx, stability_rhs(ctx), h)

    eig <- eigen(J)
    # Sort by decreasing real part (most unstable first).
    ord   <- order(Re(eig$values), decreasing = TRUE)
    evals <- as.complex(eig$values[ord])
    evecs <- eig$vectors[, ord, drop = FALSE]

    max_real_part <- Re(evals[1])
    omega1 <- Im(evals[1])
    dominant_period <- if (abs(omega1) < 1e-10) Inf else 2 * pi / abs(omega1)

    is_complex <- abs(Im(evals)) > 1e-8
    if (any(is_complex)) {
        # The complex eigenvalue closest to (or furthest across) the imaginary
        # axis is the one whose oscillation the dynamics show. It is named after
        # what it is — the leading oscillatory mode — and not after a Hopf
        # bifurcation: whether it is anywhere near crossing the axis is a
        # question about its real part, which this function reports rather than
        # presumes. Its eigenvector is returned alongside it: it need not be
        # among the leading two, and picking it out of `leading_eigenvectors` by
        # index would silently return the shape of a different mode whenever it
        # is not.
        osc_idx <- which(is_complex)[which.max(Re(evals[is_complex]))]
        leading_oscillatory_eigenvalue  <- evals[osc_idx]
        oscillation_period      <- 2 * pi / abs(Im(leading_oscillatory_eigenvalue))
        leading_oscillatory_eigenvector <- single_eigenvector(
            stability_eigenvectors(evecs, ctx, cols = osc_idx))
    } else {
        leading_oscillatory_eigenvalue  <- NULL
        oscillation_period      <- NULL
        leading_oscillatory_eigenvector <- NULL
    }

    list(
        eigenvalues          = evals,
        max_real_part        = max_real_part,
        stable               = max_real_part < 0,
        dominant_period      = dominant_period,
        oscillation_period          = oscillation_period,
        leading_oscillatory_eigenvalue      = leading_oscillatory_eigenvalue,
        leading_oscillatory_eigenvector     = leading_oscillatory_eigenvector,
        n_active             = length(ctx$x0),
        leading_eigenvectors = stability_eigenvectors(evecs, ctx),
        params               = ctx$params
    )
}

#' Analyse the stability of mizer's numerical time step
#'
#' `r lifecycle::badge("experimental")`
#' Computes the eigenvalues \eqn{\mu_i} of the linearised one-step-ahead map at
#' the steady state stored in `params@initial_n`, for a given step size `dt`.
#' These describe how mizer's numerical scheme, rather than the model, behaves
#' near the steady state: the map does not amplify perturbations when the
#' spectral radius \eqn{\max_i|\mu_i|} is less than 1.
#'
#' This is the numerical counterpart of [getStability()], which analyses the
#' model itself and involves no time step at all. Use `getStability()` to ask
#' whether the steady state of the *model* is stable, and this function to ask
#' what mizer's solver does at a particular `dt`. The two can disagree, and that
#' disagreement is the point: the implicit transport solve damps oscillations
#' artificially, so a physically unstable steady state can have a spectral
#' radius below 1 at a large `dt`, and the simulation then sits at a state the
#' model does not actually hold.
#'
#' ## The map that is linearised
#'
#' One step is what [project()] takes with `method = "euler"`: the rates are
#' evaluated at the state at the start of the step, and the resulting transport
#' problem is solved implicitly,
#' \deqn{A(N^t, n_{pp}^t)\,N^{t+1} = S(N^t, n_{pp}^t),}
#' with the same `project_n_loop()` C++ Thomas solver and the same spatial
#' scheme ([second_order_w()]) as the regular dynamics.
#'
#' Because the rates are evaluated at the *input* state, the step is not fully
#' implicit, and the discrete eigenvalues therefore cannot be converted into
#' continuous-time eigenvalues by any exact algebraic relation. That conversion
#' is what [getStability()] avoids by differentiating the rates of change
#' themselves.
#'
#' The resource is advanced by the model's own `resource_dynamics` function, the
#' one [project()] calls. Nothing is substituted for it: the map that is
#' differentiated here reproduces a single `project(method = "euler")` step
#' exactly, which is what makes the spectral radius a statement about mizer's
#' solver rather than about a nearby scheme.
#'
#' @inheritParams getStability
#' @param dt The time step size of the one-step map. Default `1`.
#' @return A named list with the following components:
#'   \describe{
#'     \item{`discrete_eigenvalues`}{Complex vector of the eigenvalues
#'       \eqn{\mu_i} of the one-step map, sorted by decreasing modulus.}
#'     \item{`spectral_radius`}{\eqn{\max_i|\mu_i|}. Less than 1 means the
#'       numerical scheme is stable at this `dt`.}
#'     \item{`stable`}{Logical: `TRUE` when `spectral_radius < 1`.}
#'     \item{`dt`}{The step size the map was evaluated at.}
#'     \item{`n_active`}{Dimension of the Jacobian.}
#'     \item{`leading_eigenvectors`}{The eigenvectors of the two
#'       largest-modulus eigenvalues, in the same shape as for
#'       [getStability()].}
#'     \item{`params`}{The validated `params` object the analysis was made at.}
#'   }
#' @inheritSection getStability Requires smooth dynamics
#' @seealso [getStability()], [findSteadyState()]
#' @export
getDiscreteStability <- function(params,
                                 effort = params@initial_effort,
                                 h = 1e-4,
                                 dt = 1) {
    assert_that(is.number(h), is.finite(h), h > 0,
                is.number(dt), is.finite(dt), dt > 0)
    ctx <- stability_context(params,
                             effort = effort,
                             fn_name = "getDiscreteStability",
                             map_name = "the one-step map")
    L <- stability_jacobian(ctx, stability_step(ctx, dt), h)

    eig <- eigen(L)
    ord   <- order(Mod(eig$values), decreasing = TRUE)
    evals <- as.complex(eig$values[ord])
    evecs <- eig$vectors[, ord, drop = FALSE]
    spectral_radius <- Mod(evals[1])

    list(
        discrete_eigenvalues = evals,
        spectral_radius      = spectral_radius,
        stable               = spectral_radius < 1,
        dt                   = dt,
        n_active             = length(ctx$x0),
        leading_eigenvectors = stability_eigenvectors(evecs, ctx),
        params               = ctx$params
    )
}
