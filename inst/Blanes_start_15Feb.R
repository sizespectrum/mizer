species_params <- read.csv("vignettes/Blanes_species_params_full.csv", sep = ";")
library(readr)
theta <- read_csv("vignettes/Blanes_theta_full3.csv", 
                  col_names = FALSE)
theta <- t(as.matrix(theta))

## Choose B ----

B <- c(detritus = 1, carrion = 1)

## Choose rho ----
my_rho  <- 10^20
no_sp <- dim(species_params)[1]
rho <- array(0, dim = c(no_sp,2))
#rho[1,2] <- my_rho
rho[,1] <- species_params$detQ*10^20
rho[,2] <- species_params$scavQ*10^20
################################

## Set up params ----
# with constant resource biomass
resource_dynamics <- 
    list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"],
         "carrion" = function(params, n, n_pp, B, rates, dt, ...) B["carrion"])
resource_params <- list("detritus_external" = 0,
                          "detritus_proportion" = 0.5, 
                          "carrion_external" = 0)

params <- MizerParams(species_params, rho = rho,
                      resource_dynamics = resource_dynamics,
                      resource_params = resource_params)

# TODO: choose the initial_n to correspond somewhat to the actual abundances
# observed in the ecosystem.

# and constant reproduction
rdd <- getRDD(params)
params@srr <- function(rdi, species_params) {rdd}

## Run to steady state ----
sim <- project(params, B = B, t_max = 100)
plotBiomass(sim)

no_t <- dim(sim@n)[1]
n <- sim@n[no_t, , ]
n_pp <- sim@n_pp[no_t, ]

## compute carrion consumption for steady state ----
r <- getRates(params, n, n_pp, B = B)
carrion_cons <-
    B["carrion"] * sum((params@rho[, "carrion", ] * n * (1 - r$feeding_level)) %*%
                params@dw)
detritus_cons <-
    B["detritus"] - detritus_dynamics(params, n, n_pp, B, rates = r, dt = 0.1)

## Update params object
p <- params
p@resource_dynamics <- list("detritus" = detritus_dynamics,
                            "carrion" = carrion_dynamics)
p@resource_params <- list("detritus_external" = detritus_cons,
                          "detritus_proportion" = 0, 
                        "carrion_external" = carrion_cons)
p@initial_n <- n
p@initial_n_pp <- n_pp

# Retune the values of erepro so that we get the correct level of
# recruitment
mumu <- getMort(p, p@initial_n, p@initial_n_pp, B = B, effort = 0)
gg <- getEGrowth(p, p@initial_n, p@initial_n_pp, B = B)
rdd <- getRDD(p, p@initial_n, p@initial_n_pp, B = B)
# TODO: vectorise this
for (i in (1:no_sp)) {
    gg0 <- gg[i, p@w_min_idx[i]]
    mumu0 <- mumu[i, p@w_min_idx[i]]
    DW <- p@dw[p@w_min_idx[i]]
    p@species_params$erepro[i] <- p@species_params$erepro[i] *
        (p@initial_n[i, p@w_min_idx[i]] *
             (gg0 + DW * mumu0)) / rdd[i]
}

sim2 <- project(p, dt = 0.1, t_max = 100, initial_B = B)
plotBiomass(sim2)
plot(sim2)

plot(sim2@B[, "carrion"], type = "l")



###### got a steady state ##############

# allow many carrior eaters

# allow detritus eaters

# turn on detritus_proportion, and steady that
