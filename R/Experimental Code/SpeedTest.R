
source('./R/MizerParams-classNewSlots.R.R')
#source('./R/MizerParams-class.R')
source('./R/MizerSim-class.R')
#source('./R/MizerSim-classNewSlots.R')
#source('./R/project_methodsFFT.R')
 source('./R/project_methods.R')
source('./R/selectivity_funcs.R')
source('./R/summary_methods.R')
source('./R/wrapper_functions.R')
source('./R/plots.R')
source('./R/project.R')
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)

params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
params <- MizerParams(params_data, interaction = inter, no_w = 200)


#Read mizer project code
source('./R/project_methods.R')
#start clock
ptm <- proc.time() 
# run simulation
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
# Stop the clock
timeData<- proc.time() - ptm
# Print elapsed time
mizerTime <- timeData[[3]]

#Read fft project code
source('./R/project_methodsFFT2.R')
#start clock
ptm <- proc.time() 
# run simulation
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
# Stop the clock
timeData<- proc.time() - ptm
# Print elapsed time
fftTime <- timeData[[3]]

c(mizerTime, fftTime)

