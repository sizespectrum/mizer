params_data <- read.csv("vignettes/Blanes_species_params_play.csv", sep = ";")

## Choose B ----

B <- c(0, 1)
names(B) <- c("detritus", "carrion")

## Choose rho ----
my_rho  <- 10^20
no_sp <- dim(params_data)[1]
rho.array <- array(0, dim = c(no_sp,2))
rho.array[1,2] <- my_rho

## Set up params ----
# with constant resource biomass
identity_fun <- function(params, n, n_pp, B, rates, dt, ...) {return(B)}

resource_params <- list("detritus_external" = 0, 
                        "detritus_proportion" = 0, 
                        "carrion_external" = 0)
params <- multispeciesParams(params_data, rho = rho.array,
                             resource_dyn = identity_fun, 
                             resource_params = resource_params, 
                             resource_names = c("detritus","carrion"))

# TODO: choose the initial_n to correspond somewhat to the actual abundances
# observed in the ecosystem.

# and constant reproduction
rdd <- getRDD(params, params@initial_n, params@initial_n_pp, B)
params@srr <- function(rdi, species_params) {rdd}

## Run to steady state ----
sim <- project(params, B = B, t_max = 100)
plotBiomass(sim)

no_t <- dim(sim@n)[1]
n <- sim@n[no_t, , ]
n_pp <- sim@n_pp[no_t, ]

## compute carrion consumption for steady state ----
feeding_level <- getFeedingLevel(params, n, n_pp, B = B)
carrion_cons <-
    B[2] * sum((params@rho[, "carrion", ] * n * (1 - feeding_level)) %*%
                params@dw)

## Update params object
p <- params
p@resource_dyn <- dead_matter_dyn
p@resource_params <- list("detritus_external" = 0, "detritus_proportion" = 0, 
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

plot(sim2@B[, 2], type = "l")



###### got a steady state ##############

# allow many carrior eaters

# allow detritus eaters

# turn on detritus_proportion, and steady that
