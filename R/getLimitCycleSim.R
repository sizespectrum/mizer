# Limit-cycle simulation from linear stability analysis ----------------------
#
# `getLimitCycleSim()` takes the output of `findSteadyState()` (with stability
# analysis attached) and constructs a `MizerSim` covering one period of the
# limit cycle in the *linear approximation*.  The approximation is exact at a
# Hopf bifurcation (Re(lambda) = 0) and remains a good first picture of the
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
#' Near a Hopf bifurcation the dynamics are dominated by a complex-conjugate
#' pair of eigenvalues \eqn{\lambda = \sigma \pm i\omega} of the Jacobian, with
#' \eqn{\sigma} small and period \eqn{T = 2\pi/|\omega|}. `getStability()`
#' returns that pair as `hopf_eigenvalue` and its eigenvector as
#' `hopf_eigenvector`. The linearised perturbation of the full state
#' \eqn{x = (N, n_{pp})} is
#' \deqn{\delta x(t) = A\,\operatorname{Re}[e^{i\omega t}\,\mathbf{v}],}
#' where \eqn{\mathbf{v}} is that eigenvector and \eqn{A} is chosen so that the
#' largest *relative swing in species biomass* equals `amplitude`. Biomass is a
#' linear functional of the abundance, so
#' \deqn{B_i(t) = B_i^* + A\,\operatorname{Re}[e^{i\omega t} c_i], \qquad
#'       c_i = \int v_i(w)\, w \, dw,}
#' and species \eqn{i} departs from its steady biomass by at most
#' \eqn{A|c_i|}. \eqn{A} is set so that \eqn{\max_i A|c_i|/B_i^*} is
#' `amplitude`: no species' biomass moves further than that fraction from its
#' steady value, and the one that oscillates hardest moves exactly that far.
#' The integral uses [sizeIntegral()], so it follows the model's own quadrature
#' scheme and agrees with [getBiomass()] bin for bin.
#'
#' The cap is on the species that swings hardest rather than on the community
#' total, because species oscillating out of phase cancel in the total: a
#' modest total swing can be produced by wild swings in the individual species.
#'
#' The state at each time is
#' \deqn{x(t) = \max(x^* + \delta x(t),\; 0).}
#' Because the cap is on biomass, an individual size class can still be driven
#' negative while the biomass it belongs to moves only a little — a cohort
#' trough is a much larger relative excursion than the biomass integral over
#' it. That clipping is reported when it happens, and means the picture is no
#' longer the linear mode.
#'
#' The fish and resource blocks of \eqn{\mathbf{v}} carry a single common
#' normalisation, so the same \eqn{A} drives both and the resource oscillates
#' with the amplitude and phase *the mode gives it* — generally neither in step
#' with the fish nor slaved to them. `amplitude` is set on the fish biomass, so
#' how far the resource moves is a property of the mode rather than something
#' you choose.
#'
#' The growth of the mode is deliberately dropped: \eqn{e^{\sigma t}} is
#' omitted so that the cycle closes after one period. That is exact only at the
#' bifurcation itself, where \eqn{\sigma = 0}; away from it the returned object
#' shows the *shape* of the oscillation, not its envelope. \eqn{\sigma} is
#' recorded in the result's `sim_params` as `growth_rate`.
#'
#' The returned \linkS4class{MizerSim} has times running from 0 to \eqn{T}
#' (the period in time steps, typically years).
#'
#' @param x A \linkS4class{MizerParams} object at a steady state,
#'   typically the output of [findSteadyState()], or the list returned by
#'   [getStability()].
#' @param amplitude Largest relative swing in species biomass across the cycle,
#'   \eqn{\max_i \max_t |B_i(t) - B_i^*| / B_i^*}.  Default `0.1`, meaning the
#'   most strongly oscillating species departs 10\% from its steady biomass.
#' @param t_save The time interval between saved time steps in the returned
#'   \linkS4class{MizerSim}.  Defaults to `0.1`.
#' @param ... Additional arguments forwarded to [getStability()] when `x` is a
#'   `MizerParams` object.
#' @return A \linkS4class{MizerSim} object whose time axis spans one period
#'   \eqn{[0, T]} of the linearised limit cycle.
#' @seealso [getStability()], [findSteadyState()]
#' @export
getLimitCycleSim <- function(x, amplitude = 0.1, t_save = 0.1, ...) {
    # ------------------------------------------------------------------
    # 1. Get (or compute) stability analysis
    # ------------------------------------------------------------------
    if (is(x, "MizerParams")) {
        params <- x
        # `getStability()` raises the not-at-steady-state warning itself.
        stab <- getStability(params, ...)
    } else if (is.list(x) && !is.null(x$eigenvalues) && is(x$params, "MizerParams")) {
        stab <- x
        params <- stab$params
    } else {
        stop("The first argument must be a MizerParams object or the stability list returned by getStability().")
    }

    if (is.null(stab$hopf_eigenvalue)) {
        stop("No oscillatory mode detected: all eigenvalues are real. ",
             "A Hopf limit cycle requires at least one complex eigenvalue ",
             "pair. Try a parameter setting closer to a Hopf bifurcation.")
    }

    # ------------------------------------------------------------------
    # 2. The eigenvalue / eigenvector to use
    #
    # `getStability()` picks out the dominant *oscillatory* mode and returns
    # its eigenvector alongside it. That mode need not be among the leading
    # eigenvectors at all, so it has to travel with its own vector rather than
    # be looked up by index.
    # ------------------------------------------------------------------
    lam1   <- stab$hopf_eigenvalue
    v_fish <- stab$hopf_eigenvector$fish
    v_npp  <- stab$hopf_eigenvector$resource

    theta    <- Im(lam1)               # angular frequency, per year
    T_period <- 2 * pi / abs(theta)    # period in years

    # ------------------------------------------------------------------
    # 3. Time grid
    # ------------------------------------------------------------------
    n_steps <- ceiling(T_period / t_save)
    t_seq <- seq(0, by = t_save, length.out = n_steps + 1)

    # ------------------------------------------------------------------
    # 4. Amplitude scaling
    #    A * max_w( |v(w)| / N*(w) ) = amplitude  =>  A = amplitude / that max
    # ------------------------------------------------------------------
    N_ss       <- params@initial_n
    npp_ss     <- params@initial_n_pp
    n_other_ss <- params@initial_n_other

    # `amplitude` is set on species biomass, which is the quantity a reader of
    # plotBiomass() actually sees. Biomass is a linear functional of the
    # abundance, so the biomass of the perturbation is the same integral
    # applied to the eigenvector:
    #   B_i(t) = B_i* + A Re[e^{i omega t} c_i],   c_i = int v_i(w) w dw,
    # and the largest deviation species i reaches over the cycle is A |c_i|.
    # The integral goes through sizeIntegral(), so it follows the model's own
    # quadrature scheme and matches getBiomass() bin for bin. It is real and
    # linear, so the real and imaginary parts integrate separately.
    biomass_of <- function(n) {
        as.numeric(sizeIntegral(params, weighting = params@w, n = n))
    }
    B_ss  <- biomass_of(N_ss)
    c_sp  <- complex(real      = biomass_of(Re(v_fish)),
                     imaginary = biomass_of(Im(v_fish)))

    ok <- B_ss > 0 & is.finite(B_ss)
    if (!any(ok) || max(Mod(c_sp[ok]) / B_ss[ok]) <= 0) {
        stop("The leading oscillatory mode moves no species' biomass, so ",
             "there is no cycle to draw at a given biomass amplitude.")
    }
    # Cap the species that swings hardest, rather than the community total:
    # species oscillating out of phase cancel in the total, which would let a
    # modest total swing be produced by wild swings in the individual species.
    A_scale <- amplitude / max(Mod(c_sp[ok]) / B_ss[ok])

    # ------------------------------------------------------------------
    # 5. Build the MizerSim
    # ------------------------------------------------------------------
    sim <- MizerSim(params, t_dimnames = t_seq)
    sim@sim_params <- list(
        method      = "limit_cycle_linear_approx",
        period      = T_period,
        amplitude   = amplitude,
        eigenvalue  = lam1,
        growth_rate = Re(lam1)
    )

    has_other <- length(params@other_dynamics) > 0
    clipped   <- FALSE

    for (t_idx in seq_along(t_seq)) {
        t <- t_seq[t_idx]
        phase <- exp(1i * theta * t)

        # Fish: N(t) = N* + A * Re[ e^{iωt} * v_fish ]
        N_t <- N_ss + A_scale * Re(phase * v_fish)
        # Resource: the same A and the same phase factor, applied to the
        # resource block of the same eigenvector.
        npp_t <- npp_ss + A_scale * Re(phase * v_npp)

        # A cell that is zero in the steady state can be pushed below zero by
        # any perturbation at all, and clipping it is routine rather than a
        # sign that the amplitude is too large. Only positive cells count.
        # A tolerance, because `amplitude = 1` puts the extreme cell exactly
        # at zero and rounding decides the sign.
        clipped <- clipped ||
            any(N_t[N_ss > 0] < -1e-10 * N_ss[N_ss > 0]) ||
            any(npp_t[npp_ss > 0] < -1e-10 * npp_ss[npp_ss > 0])
        sim@n[t_idx, , ]  <- pmax(N_t, 0)
        sim@n_pp[t_idx, ] <- pmax(npp_t, 0)

        # Effort: constant at initial effort
        sim@effort[t_idx, ] <- params@initial_effort

        # Other components: constant at steady state
        if (has_other) {
            sim@n_other[t_idx, ] <- n_other_ss
        }
    }

    if (clipped) {
        signal_info("amplitude", paste0(
            "An amplitude of ", signif(amplitude, 3), " drives some ",
            "abundances negative, and they have been clipped at zero, so the ",
            "cycle shown is no longer the linear mode. Reduce `amplitude` to ",
            "stay inside the linear approximation."),
            level = 1, severity = "warning", unhandled = "show",
            class = "info_limit_cycle_clipped")
    }

    sim
}
