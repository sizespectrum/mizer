library(mizer)

data(NS_species_params_gears)
data(inter)

my.array <- array((1:24), dim=c(12,2))


resource_params <- list("detritus_external" = .1, "detritus_proportion" = 0.1, "carrion_external"=1)

params <- multispeciesParams(NS_species_params_gears, inter, rho = my.array, 
                             resource_dyn = dead_matter_dyn, resource_params = resource_params, 
                             resource_names= c("detritus","carrion"))

sim <- project(params, dt=0.1, t_save = 1, t_max = 100)

plot(sim@B[,1])

plot(sim@B[,2])


#params	
#n	
#n_pp	
#B	
#dt	

#resource_params	
#Either NULL (default) or a list of parameters needed by the resource_dyn function

#params@resource_params$carrion_external.

#params@resource_params$detritus_external
#params@resource_params$detritus_proportion
