sp <- read.csv("inst/blanes/species_params.csv", sep = ";")
sp$k_vb[1:2] <- sp$k_vb[1:2] * 5
sp$k_vb[c(8, 10)] <- sp$k_vb[c(8, 10)] * 4

sp$interaction_p = rep(0, no_sp),

theta <- read.csv("inst/blanes/theta.csv", header = FALSE)
theta <- t(as.matrix(theta))

B <- c(detritus = 1, carrion = 1)

q <- 0.8
n <- 0.7
p <- 0.7
lambda <- 2 + q - n
kappa <- 1

sp$alpha <- 0.1
# default for h and ks
f0 <- 0.6
sp$h <- 3 * sp$k_vb / (sp$alpha * f0) * sp$w_inf^(1/3)
sp$ks <- 0 # 0.2 * sp$h * sp$alpha * f0
# default for rho
no_sp <- nrow(sp)
rho <- array(0, dim = c(no_sp, 2))
rho[, 1] <- sp$h * f0 / (1 - f0) * sp$detQ
rho[, 2] <- sp$h * f0 / (1 - f0) * sp$scavQ
# default for gamma
lm2 <- lambda - 2
ae <- sqrt(2 * pi) * sp$sigma * sp$beta^lm2 *
    exp(lm2^2 * sp$sigma^2 / 2) *
    # The factor on the following lines takes into account the cutoff
    # of the integral at 0 and at beta + 3 sigma
    (pnorm(3 - lm2 * sp$sigma) + 
         pnorm(log(sp$beta)/sp$sigma + 
                   lm2 * sp$sigma) - 1)
sp$gamma <- sp$h / (kappa * ae) * f0 / (1 - f0)

# Set up params with constant resource biomass
resource_dynamics <- 
    list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"],
         "carrion" = function(params, n, n_pp, B, rates, dt, ...) B["carrion"])

params <- MizerParams(sp, rho = rho,
                      interaction = theta,
                      resource_dynamics = resource_dynamics,
                      q = q, p = p, n = n, lambda = lambda, kappa = kappa,
                      f0 = f0, z0pre = 2)
params@initial_B <- B
params@initial_n <- params@initial_n * 100
params@initial_n_pp[] <- 0
params@cc_pp[] <- 0

## Construct initial size spectra ----

# Get constants for steady-state solution
mu0 <- (1 - f0) * sqrt(2 * pi) * kappa * sp$gamma * sp$sigma *
    (sp$beta ^ (n - 1)) * exp(sp$sigma ^ 2 * (n - 1) ^ 2 / 2)
hbar <- sp$alpha * sp$h * f0 - sp$ks
if (any(hbar < 0)) {
    stop("The feeding level is not sufficient to maintain the fish.")
}
pow <- mu0 / hbar / (1 - n)
if (any(pow < 1)) {
    message("The ratio of death rate to growth rate is too small, leading to
                an accumulation of fish at their largest size.")
}

w <- params@w
dw <- params@dw
for (i in 1:no_sp) {
    mumu <- mu0[i] * w^(n - 1)  # Death rate
    gg <- hbar[i] * w^n * (1 - params@psi[i, ])  # Growth rate
    w_inf_idx <- sum(w <= sp$w_inf[i])
    idx <- params@w_min_idx[i]:(w_inf_idx - 1)
    idxs <- params@w_min_idx[i]:(w_inf_idx)
    # Steady state solution of the upwind-difference scheme used in project
    params@initial_n[i, idxs] <- params@initial_n[i, params@w_min_idx[i] + 1] * 
        c(1, cumprod(gg[idx] / ((gg + mumu * dw)[idx + 1])))
}

#plotlySpectra(params, plankton = FALSE, total = TRUE)

# params@sc <- kappa * w^(-lambda)
# params <- retuneAbundance(params)
# g <- getEGrowth(params)[23,]
# gr <- getEReproAndGrowth(params)[23, ]
# r <- getERepro(params)[23,]
# plot(params@w, gr, type = "l", log = "x")
# lines(params@w, g, col = "blue")
# lines(params@w, r, col = "red")

# and constant reproduction
rdd <- getRDD(params)
params@srr <- function(rdi, species_params) {rdd}

## Run to steady state ----
sim <- project(params, t_max = 200)
#plotlyBiomass(sim)
plotlySpectra(sim, plankton = FALSE, total = TRUE)

no_t <- dim(sim@n)[1]
n <- sim@n[no_t, , ]
n_pp <- sim@n_pp[no_t, ]

## compute carrion consumption for steady state ----
r <- getRates(params, n, n_pp, B)
(r$rdi - r$rdd) / r$rdi

# carrion_cons <-
#     B["carrion"] * sum((params@rho[, "carrion", ] * n * (1 - r$feeding_level)) %*%
#                 params@dw)
# detritus_cons <-
#     B["detritus"] - detritus_dynamics(params, n, n_pp, B, rates = r, dt = 0.1)

## Update params object
p <- params
# p@resource_dynamics <- list("detritus" = detritus_dynamics,
#                             "carrion" = carrion_dynamics)
p@resource_params <- list("detritus_external" = 0, #detritus_cons,
                          "detritus_proportion" = 0,
                        "carrion_external" = 0) #carrion_cons)
p@initial_n <- n
p@initial_n_pp <- n_pp

p@srr <- srrBevertonHolt

# Retune the values of erepro so that we get the correct level of
# recruitment
p@species_params$erepro <- r$rdd / r$rdi

saveRDS(p, file = "inst/tuning/params.rds")

# Recompute all species
mumu <- getMort(p, effort = 0)
gg <- getEGrowth(p)
for (sp in 1:length(p@species_params$species)) {
    w_inf_idx <- sum(p@w < p@species_params[sp, "w_inf"])
    idx <- p@w_min_idx[sp]:(w_inf_idx - 1)
    if (any(gg[sp, idx] == 0)) {
        stop("Can not compute steady state due to zero growth rates")
    }
    n0 <- p@initial_n[sp, p@w_min_idx[sp]]
    p@initial_n[sp, ] <- 0
    p@initial_n[sp, p@w_min_idx[sp]:w_inf_idx] <- 
        c(1, cumprod(gg[sp, idx] / ((gg[sp, ] + mumu[sp, ] * p@dw)[idx + 1]))) *
        n0
}

plotlySpectra(p, plankton = FALSE, total = TRUE)

# Retune the values of erepro so that we get the correct level of
# recruitment
mumu <- getMort(p, effort = 0)
gg <- getEGrowth(p)
rdd <- getRDD(p)
# TODO: vectorise this
for (i in (1:length(p@species_params$species))) {
    gg0 <- gg[i, p@w_min_idx[i]]
    mumu0 <- mumu[i, p@w_min_idx[i]]
    DW <- p@dw[p@w_min_idx[i]]
    p@species_params$erepro[i] <- p@species_params$erepro[i] *
        p@initial_n[i, p@w_min_idx[i]] *
        (gg0 + DW * mumu0) / rdd[i]
}

sim2 <- project(p, dt = 0.1, t_max = 100)
plotlyBiomass(sim2)
plot(sim2)

plot(sim2@B[, "carrion"], type = "l")


