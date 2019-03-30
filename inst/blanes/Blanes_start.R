params_data <- read.csv("C:/Users/user/Desktop/mizer projects/detritus/mizer/vignettes/Blanes_species_params_play.csv", sep = ";")

params <- multispeciesParams(params_data)

sim <- project(params, dt=0.1, t_save = .1, t_max = 100)
plot(sim)
