# Limit-cycle simulation from linear stability analysis ----------------------
#
# `getLimitCycleSim()` takes the output of `steadyNewton()` (with stability
# analysis attached) and constructs a `MizerSim` covering one period of the
# limit cycle in the *linear approximation*.  The approximation is exact at a
# Hopf bifurcation (|lambda| = 1) and remains a good first picture of the
# oscillation pattern for parameter values close to the bifurcation.

#' Construct a MizerSim of the linearised limit cycle
#'
#' `r lifecycle::badge("experimental")`
#' Using the leading complex eigenvector from [getStability()], constructs a
#' \linkS4class{MizerSim} object covering one period of the limit cycle in the
#' linear approximation.  The result can be inspected with all standard mizer
#' plotting functions (e.g. [plotBiomass()], [plotSpectra()]).
#'
#' ## Mathematical background
#'
#' Near a Hopf bifurcation the dominant eigenvalue of the linearised one-step
#' map is a complex conjugate pair \eqn{\lambda = r e^{i\theta}} with
#' \eqn{r \approx 1} and angular frequency \eqn{\theta = 2\pi/T} per time
#' step. The linearised perturbation is
#' \deqn{\delta N(t) = A\,\operatorname{Re}[e^{i\theta t}\,\mathbf{v}],}
#' where \eqn{\mathbf{v}} is the leading complex eigenvector (normalised so
#' \eqn{\max_w |\mathbf{v}(w)| = 1}) and \eqn{A} is chosen so that the
#' maximum *relative* perturbation
#' \eqn{\max_{t,w}|\delta N(t,w)|/N^*(w) = } `amplitude`.
#' The full state at each time step is
#' \deqn{N(t) = \max(N^* + \delta N(t),\; 0)}
#' (clipping prevents negative abundances at large amplitudes).
#'
#' The resource is placed at its quasi-static semichemostat equilibrium
#' \eqn{n_{pp}^*(N(t))} at each step, consistent with the reduced Jacobian
#' used by [getStability()].
#'
#' The returned \linkS4class{MizerSim} has times running from 0 to \eqn{T}
#' (the period in time steps, typically years).
#'
#' @param params A \linkS4class{MizerParams} object at a steady state,
#'   typically the output of [steadyNewton()].  If `attr(params, "stability")`
#'   exists and contains `leading_eigenvectors` it is used directly; otherwise
#'   [getStability()] is called.
#' @param amplitude Maximum relative perturbation
#'   \eqn{\max_w |\delta N(t,w)|/N^*(w)} across the limit cycle.  Default
#'   `0.1` (10\%).
#' @param n_points Number of time points in the returned
#'   \linkS4class{MizerSim}.  Defaults to `ceiling(T) + 1`, giving one point
#'   per time step (year) over exactly one period.  Increase for smoother
#'   plots, e.g. `n_points = 200`.
#' @param ... Additional arguments forwarded to [getStability()] when
#'   stability has not already been computed and cached on `params`.
#' @return A \linkS4class{MizerSim} object whose time axis spans one period
#'   \eqn{[0, T]} of the linearised limit cycle.
#' @seealso [getStability()], [steadyNewton()]
#' @export
getLimitCycleSim <- function(params, amplitude = 0.1, n_points = NULL, ...) {
    # ------------------------------------------------------------------
    # 1. Get (or compute) stability analysis
    # ------------------------------------------------------------------
    stab <- attr(params, "stability")
    if (is.null(stab) || is.null(stab$leading_eigenvectors)) {
        message("Computing stability analysis ...")
        stab <- getStability(params, ...)
    }

    if (is.null(stab$hopf_period)) {
        stop("No oscillatory mode detected: all eigenvalues are real. ",
             "A Hopf limit cycle requires at least one complex eigenvalue pair.")
    }

    # ------------------------------------------------------------------
    # 2. Identify the eigenvalue / eigenvector to use
    # ------------------------------------------------------------------
    lam1 <- stab$eigenvalues[1]

    if (abs(Im(lam1)) <= 1e-8) {
        # Dominant mode is real (monotone); the oscillatory Hopf mode is
        # not the leading one.  Error with guidance.
        stop("The dominant eigenvalue is real (monotone dynamics). ",
             "The oscillatory Hopf mode (period ",
             round(stab$hopf_period, 1), " time steps) is not the leading ",
             "eigenvalue; getLimitCycleSim() requires the Hopf mode to be ",
             "dominant.  Try using a parameter setting closer to the Hopf ",
             "bifurcation, where the complex pair has larger modulus.")
    }

    # Leading eigenvectors are stored as (sp × w × 2) array, or as a list
    # with $fish when include_resource = TRUE was used.
    lev <- stab$leading_eigenvectors
    v_use <- if (is.array(lev)) lev[, , 1] else lev$fish[, , 1]

    theta    <- Arg(lam1)               # angular frequency per time step
    T_period <- 2 * pi / abs(theta)    # period in time steps

    # ------------------------------------------------------------------
    # 3. Time grid
    # ------------------------------------------------------------------
    if (is.null(n_points)) {
        t_seq <- seq(0, ceiling(T_period))   # integer steps, one per year
    } else {
        t_seq <- seq(0, T_period, length.out = as.integer(n_points))
    }

    # ------------------------------------------------------------------
    # 4. Amplitude scaling
    #    A * max_w( |v(w)| / N*(w) ) = amplitude  =>  A = amplitude / that max
    # ------------------------------------------------------------------
    N_ss       <- params@initial_n
    npp_ss     <- params@initial_n_pp
    n_other_ss <- params@initial_n_other

    active  <- N_ss > 0 & is.finite(N_ss)
    rel_mod <- ifelse(active,
                      Mod(v_use) / pmax(N_ss, .Machine$double.eps),
                      0)
    A_scale <- amplitude / max(rel_mod)

    # ------------------------------------------------------------------
    # 5. Build the MizerSim
    # ------------------------------------------------------------------
    sim <- MizerSim(params, t_dimnames = t_seq)
    sim@sim_params <- list(
        method     = "limit_cycle_linear_approx",
        period     = T_period,
        amplitude  = amplitude,
        eigenvalue = lam1
    )

    has_other <- length(params@other_dynamics) > 0

    for (t_idx in seq_along(t_seq)) {
        t <- t_seq[t_idx]

        # Fish: N(t) = N* + A * Re[ e^{iθt} * v ]
        delta_N <- A_scale * Re(exp(1i * theta * t) * v_use)
        N_t     <- pmax(N_ss + delta_N, 0)
        sim@n[t_idx, , ] <- N_t

        # Resource: quasi-static semichemostat equilibrium at N_t
        if (params@resource_dynamics == "resource_semichemostat") {
            sim@n_pp[t_idx, ] <- resource_steady_semichemostat(params, N_t,
                                                               n_other_ss)
        } else {
            sim@n_pp[t_idx, ] <- npp_ss
        }

        # Effort: constant at initial effort
        sim@effort[t_idx, ] <- params@initial_effort

        # Other components: constant at steady state
        if (has_other) {
            sim@n_other[t_idx, ] <- n_other_ss
        }
    }

    sim
}
