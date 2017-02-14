source('./R/MizerParams-class.R')
source('./R/MizerSim-class.R')
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

#params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7,
#                              sigma = 2, no_w=100)

#params <- set_trait_model(no_sp = 2, min_w_inf = 10, max_w_inf = 1e5)
#sim <- project(params)

params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
params <- MizerParams(params_data, interaction = inter)
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
plot(sim)
object <- sim@params 

nt <- dim(sim@n_pp)[1]  # index of last time step
n_pp <- sim@n_pp[nt, ]
n <- sim@n[nt, , ]


gPP <- getPhiPrey(object, n, n_pp)
i <- 6
plot(log(object@w), gPP[i,])

##################################

source('./R/MizerParams-class.R')
source('./R/MizerSim-class.R')
source('./R/project_methodsFFT.R')
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

#params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7,
#                              sigma = 2, no_w=100)

#params <- set_trait_model(no_sp = 2, min_w_inf = 10, max_w_inf = 1e5)
#sim <- project(params)

params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
params <- MizerParams(params_data, interaction = inter)
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
#plot(sim)
object <- sim@params 

nt <- dim(sim@n_pp)[1]  # index of last time step
n_pp <- sim@n_pp[nt, ]
n <- sim@n[nt, , ]


gPP <- getPhiPrey(object, n, n_pp)
i <- 6
lines(log(object@w), gPP[i,])