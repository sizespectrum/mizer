# Our first attempt at setting up a mizer model for the lagoon
species_params <- read.csv("inst/algarve/parameters for mizer.csv")
species_params$r_max <- Inf
interaction <- read.csv("inst/algarve/interation matrix CSV.csv", row.names = 1)
interaction <- as.matrix(interaction)
interaction_p <- rep(0, nrow(species_params))
names(interaction_p) <- species_params$species
interaction_p["Engraulis encrasicholus"] <- 1
interaction_p["Sardina pilchardus"] <- 1
interaction_p

species_params$h <- 10^20
species_params$ks <- 0

rho <- rep(0, nrow(species_params))
rho[20] <- 1
rho <- array(rho, dim = c(20, 1))

resource_dynamics <- list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"])

params <- MizerParams(species_params,
                      interaction = interaction,
                      interaction_p = interaction_p,
                      rho = rho,
                      resource_dynamics = resource_dynamics)
params@initial_B <- 10^8
names(params@initial_B) <- "detritus"
plotSpectra(params, plankton = FALSE)

# Run to steady state with constant reproduction
rdd <- getRDD(params)
params@srr <- function(rdi, species_params) {rdd}
sim <- project(params, t_max = 100)
plotBiomass(sim)
