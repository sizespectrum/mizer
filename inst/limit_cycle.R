library(mizer)
library(ggplot2)

norm_box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
    phi <- rep(1, length(ppmr))
    phi[ppmr > ppmr_max | ppmr < ppmr_min] <- 0
    phi[1] <- 0
    phi / (sum(phi) * (log(ppmr[2]) - log(ppmr[1])))
}

exact_logistic_immigration_step <- function(n, rate, capacity, immigration,
                                            mortality, dt) {
    result <- n
    active <- is.finite(n) & is.finite(rate) & is.finite(capacity) &
        is.finite(immigration) & is.finite(mortality) &
        rate > 0 & capacity > 0 & immigration >= 0
    if (dt == 0 || !any(active)) return(result)
    
    n0 <- n[active]; r <- rate[active]; k <- capacity[active]
    i  <- immigration[active]; mu <- mortality[active]
    a  <- r - mu; b <- r / k
    next_n <- numeric(length(n0))
    
    # With immigration the quadratic has one positive and one negative root;
    # the alternative root formulae avoid cancellation when abs(a) is large.
    has_immigration <- i > 0
    if (any(has_immigration)) {
        idx <- which(has_immigration)
        ai <- a[idx]; bi <- b[idx]; ii <- i[idx]
        d  <- sqrt(ai^2 + 4 * bi * ii)
        n_plus  <- ifelse(ai >= 0, (ai + d) / (2 * bi), 2 * ii / (d - ai))
        n_minus <- ifelse(ai <= 0, (ai - d) / (2 * bi), -2 * ii / (ai + d))
        ratio <- ((n0[idx] - n_plus) / (n0[idx] - n_minus)) * exp(-d * dt)
        next_n[idx] <- (n_plus - ratio * n_minus) / (1 - ratio)
    }
    
    # The zero-immigration limit is handled separately, including rate == mortality.
    no_immigration <- !has_immigration
    if (any(no_immigration)) {
        idx <- which(no_immigration)
        az <- a[idx]; bz <- b[idx]; nz <- n0[idx]
        value <- numeric(length(idx))
        zero_rate <- az == 0
        value[zero_rate] <- nz[zero_rate] / (1 + bz[zero_rate] * nz[zero_rate] * dt)
        
        positive <- az > 0 & nz > 0
        phi <- -expm1(-az[positive] * dt) / az[positive]
        value[positive] <- nz[positive] / (exp(-az[positive] * dt) + bz[positive] * nz[positive] * phi)
        
        negative <- az < 0
        exp_adt <- exp(az[negative] * dt)
        phi <- expm1(az[negative] * dt) / az[negative]
        value[negative] <- nz[negative] * exp_adt / (1 + bz[negative] * nz[negative] * phi)
        next_n[idx] <- value
    }
    
    # Round-off can only create tiny negative values; the analytic solution is
    # non-negative for non-negative initial abundance and immigration.
    result[active] <- pmax(next_n, 0)
    result
}

plankton_logistic <- function(params, n, n_pp, n_other, rates, dt = 0.1, ...) {
    exact_logistic_immigration_step(
        n = n_pp, rate = params@rr_pp, capacity = params@cc_pp,
        immigration = immigration, mortality = rates$resource_mort, dt = dt
    )
}

# Build model ----
p <- list(
    dt = 0.02, method = "tr_bdf2", dx = 0.1, w_min = 0.0003, w_inf = 66.5,
    ppmr_min = 100, ppmr_max = 30000, gamma = 750, alpha = 0.85, K = 0.1,
    mu_l = 0, w_l = 0.03, rho_l = 5,
    mu_0 = 1, rho_b = -0.25,
    w_s = 0.5, rho_s = 1,
    w_mat = 10, rho_m = 15, rho_inf = 0.2, epsilon_R = 0.1,
    w_pp_cutoff = 0.1, r0 = 10, a0 = 100, i0 = 100, rho = 0.85, lambda = 2
)

theta                 <- 0.3   # fraction of baseline cannibalism strength kept
interaction_resource  <- 1     # this species' access to the background resource
knife_edge_size       <- 10    # fishing gear cutoff, in grams
kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))

species_params_df <- data.frame(
    species = "Anchovy", w_min = p$w_min, w_mat = p$w_mat, m = p$rho_inf + 2/3,
    w_inf = p$w_inf, erepro = p$epsilon_R, alpha = p$K, ks = 0, gamma = p$gamma,
    q = p$alpha, ppmr_min = p$ppmr_min, ppmr_max = p$ppmr_max,
    pred_kernel_type = "norm_box", h = Inf, R_max = Inf, linecolour = "brown",
    stringsAsFactors = FALSE
)

params <- newMultispeciesParams(
    species_params_df, no_w = round(log(p$w_inf / p$w_min) / p$dx),
    lambda = p$lambda, kappa = kappa, w_pp_cutoff = p$w_pp_cutoff#,
    #resource_dynamics = "plankton_logistic"
)
resource_rate(params) <- p$r0 * w_full(params)^(p$rho - 1)
immigration <- p$i0 * w_full(params)^(-p$lambda) * exp(-6.9 * (p$lambda - 1))

interaction_matrix(params) <- theta
species_params(params)$interaction_resource <- interaction_resource

w    <- w(params)
mu_b <- rep(0, length(w))
mu_b[w <= p$w_s] <- (p$mu_0 * (w / p$w_min)^p$rho_b)[w < p$w_s]
mu_s <- min(mu_b[w <= p$w_s])
mu_b[w >= p$w_s] <- (mu_s * (w / p$w_s)^p$rho_s)[w >= p$w_s]
mu_b <- mu_b + p$mu_l / (1 + (w / p$w_l)^p$rho_l)

mort    <- ext_mort(params)
mort[]  <- mu_b
ext_mort(params) <- mort

ss <- projectToSteady(params, t_max = 150, t_per = 0.2, dt = p$dt,
                      method = p$method, return_sim = TRUE)
plotBiomass(ss)
conv <- attr(ss, "convergence")
start <- conv$years - conv$period
pp <- getParams(ss, c(start, conv$years))

ps <- steadyNewton(pp, reproduction = "dynamic")
stab <- getStability(ps, reproduction = "dynamic", include_resource = TRUE)
lc <- getLimitCycleSim(stab)

sims <- project(ps, t_max = 10, dt = p$dt, method = p$method)
plotBiomass(sims)

effort_seq <- c(1, 2, 3, 5, 9, 20, 50, 100)
params <- finalParams(ss)

mean_yield <- sapply(effort_seq, function(effort) {
    sim <- projectToSteady(params, dt = p$dt, method = p$method, t_per = 0.2,
                           return_sim = TRUE, effort = effort)
    conv <- attr(sim, "convergence")
    start <- conv$years - conv$period
    mean(getYield(sim)[getTimes(sim) >= start, ])
})

yield_df <- data.frame(effort = effort_seq, yield = mean_yield)

ggplot(yield_df, aes(effort, yield)) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_log10() +
    labs(x = "Fishing effort (Constant, log scale)", y = "Mean yield",
         title = "Yield curve: theta=0.3, interaction_resource=1, knife_edge=10") +
    theme_minimal()
