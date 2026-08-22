library(mizer)
library(ggplot2)

norm_box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
    phi <- rep(1, length(ppmr))
    phi[ppmr > ppmr_max | ppmr < ppmr_min] <- 0
    phi[1] <- 0
    phi / (sum(phi) * (log(ppmr[2]) - log(ppmr[1])))
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

theta                 <- 0.2   # fraction of baseline cannibalism strength kept
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
stab$stable

# The eigenvalues are already the continuous-time ones, sorted by decreasing
# real part.
lambda_cont <- stab$eigenvalues
unstable_mode <- lambda_cont[1]

cat("Is continuous system stable?", stab$stable, "\n")
cat("Maximum real part:", stab$max_real_part, "\n")

if (!stab$stable && Im(unstable_mode) != 0) {
    cat("Hopf period:", stab$dominant_period, "years\n")
}

# How the numerical step at the dt used above sees the same state:
getDiscreteStability(ps, reproduction = "dynamic", include_resource = TRUE,
                     dt = p$dt)$spectral_radius

# According to this, the steady state is stable.
# Let's test that by making a small perturbation
ps_pert <- ps
initialN(ps_pert) <- initialN(ps) * 1.01
sim <- projectToSteady(ps_pert, t_per = 0.2, dt = p$dt, tol = 1e-6,
                       method = p$method, return_sim = TRUE)
plotBiomass(sim)


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
