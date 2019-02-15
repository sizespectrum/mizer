params_data <- read.csv("C:/Users/user/Desktop/mizer projects/detritus/mizer/vignettes/Blanes_species_params_play.csv", sep = ";")

params <- multispeciesParams(params_data)

sim <- project(params, dt=0.1, t_save = .1, t_max = 100)
plot(sim)

##########################


# 1 choose B to be
B <- sum(getBiomass(sim)[dim(getBiomass(sim))[1],])

# 2. get rho
Enc <- getEncounter(params,sim@n[dim(sim@n)[1],,],sim@n_pp[dim(sim@n_pp)[1],])
plot(params@w,Enc[1,],log="xy")
lines(params@w,(0.5*10^2)*params@w^params@n,col="red")

# (0.5*10^2) = rho *steady_state_biomass_of_detritus

# match makes red curve look like Enc
match <- (0.5*10^2)

my_rho  <- match/B


# # # #

#sum((params@rho[, "carrion", ] * n * (1 - rates$feeding_level)) %*%
#        params@dw)
feeding_level <- getFeedingLevel(params,sim@n[dim(sim@n)[1],,],sim@n_pp[dim(sim@n_pp)[1],])

#my_rho <- match/(steady_state_biomass_of_detritus*
#    sum((1-feeding_level[1,])*sim@n[dim(sim@n)[1],1,]*(params@w^params@n)*params@dw))

carrion_cons <- (B*
            sum((1-feeding_level[1,])*sim@n[dim(sim@n)[1],1,]*(params@w^params@n)*params@dw)) * my_rho 

# # # #


no_sp <- dim(params@species_params)[1]
rho.array <- array(0, dim=c(no_sp,2))
rho.array[1,2] <- my_rho

resource_params <- list("detritus_external" = 0, "detritus_proportion" = 0, "carrion_external" = carrion_cons)
params <- multispeciesParams(params_data, rho = rho.array, 
                             resource_dyn = dead_matter_dyn, resource_params = resource_params, 
                             resource_names= c("detritus","carrion"))
initial_B <- c(0, B)
names(initial_B) <- c("detritus","carrion")
sim <- project(params, dt=0.1, t_save = 1, t_max = 100,initial_B = initial_B)
plot(sim)

plot(sim@B[,1],type="l")

plot(sim@B[,2])

