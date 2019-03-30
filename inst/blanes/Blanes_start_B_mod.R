params_data <- read.csv("C:/Users/user/Desktop/mizer projects/detritus/mizer/vignettes/Blanes_species_params_play.csv", sep = ";")

params <- multispeciesParams(params_data)

sim <- project(params, dt=0.1, t_save = .1, t_max = 100)
plot(sim)

##########################


# 1 choose B 

B <- sum(getBiomass(sim)[dim(getBiomass(sim))[1],])

######################
# 2. get rho
Enc <- getEncounter(params,sim@n[dim(sim@n)[1],,],sim@n_pp[dim(sim@n_pp)[1],])
plot(params@w,Enc[1,],log="xy")
lines(params@w,(0.5*10^2)*params@w^params@n,col="red")

# (0.5*10^2) = rho *steady_state_biomass_of_detritus

# match makes red curve look like Enc
match <- (0.5*10^2)

my_rho  <- match/B
no_sp <- dim(params@species_params)[1]
rho.array <- array(0, dim=c(no_sp,2))
rho.array[1,2] <- my_rho

################ make intermediate object

resource_params <- list("detritus_external" = 0, "detritus_proportion" = 0, "carrion_external" = 0)
params <- multispeciesParams(params_data, rho = rho.array, 
                             resource_dyn = dead_matter_dyn, resource_params = resource_params, 
                             resource_names= c("detritus","carrion"))

params@interaction_p[1] <- 0
params@interaction[1,] <- 0



feeding_level <- getFeedingLevel(params,sim@n[dim(sim@n)[1],,],sim@n_pp[dim(sim@n_pp)[1],], B = B)

############## compute carrion cons for steady state

carrion_cons <- (B*
                     sum((1-feeding_level[1,])*sim@n[dim(sim@n)[1],1,]*(params@w^params@n)*params@dw)) * my_rho 

###### build and run

resource_params <- list("detritus_external" = 0, "detritus_proportion" = 0, "carrion_external" = carrion_cons)
params <- multispeciesParams(params_data, rho = rho.array, 
                             resource_dyn = dead_matter_dyn, resource_params = resource_params, 
                             resource_names= c("detritus","carrion"))





initial_B <- c(0, B)
names(initial_B) <- c("detritus","carrion")


params@interaction_p[1] <- 0
params@interaction[1,] <- 0


sim <- project(params, dt=0.1, t_save = 1, t_max = 1000,initial_B = initial_B)
plot(sim)

plot(sim@B[,1],type="l")

plot(sim@B[,2], ylim=c(0,max(sim@B[,2])))

head(sim@B[,2])

###### got a steady state ##############

# allow many carrior eaters

# allow detritus eaters

# turn on detritus_proportion, and steady that
