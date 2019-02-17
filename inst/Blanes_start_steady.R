species_params <- read.csv("vignettes/Blanes_species_params_play.csv", sep = ";")

## Choose B ----

B <- c(detritus = 0, carrion = 1)

## Choose rho ----
my_rho  <- 10^20
no_sp <- dim(species_params)[1]
rho <- array(0, dim = c(no_sp,2))
rho[1,2] <- my_rho

resource_dynamics <- list(detritus = detritus_dynamics,
                            carrion = carrion_dynamics)
resource_params <- list("detritus_external" = 0, "detritus_proportion" = 0, 
                          "carrion_external" = 0)

params <- multispeciesParams(species_params, rho = rho,
                             resource_dynamics = resource_dynamics,
                             resource_params = resource_params)
params <- steady(params, t_max = 120)

sim <- project(params)
plotBiomass(sim)
