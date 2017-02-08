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

params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7,
                              sigma = 2, no_w=100)
sim <- project(params)
object <- sim@params 

# extract n_pp and n from sim object 
nt <- dim(sim@n_pp)[1]  # index of last time step
n_pp <- sim@n_pp[nt, ]
n <- sim@n[nt, , ]

#(1-feeding_level)*object@search_vol*n_total_in_size_bins

#     from line 213 of project_methods.R
#   pred_rate <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)

#    from line 297 of project_methods.R
#     m2 <- t(object@interaction) %*% colSums(aperm(pred_rate, c(2,1,3)),dims=1)[,idx_sp]


############## more about getm2 is on lines 286-300, are there are help comments above

