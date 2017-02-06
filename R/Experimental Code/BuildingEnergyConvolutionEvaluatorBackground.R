source('./R/MizerParams-class.r')
source('./R/MizerSim-class.r')
source('./R/project_methods.r')
source('./R/selectivity_funcs.r')
source('./R/summary_methods.r')
source('./R/wrapper_functions.R')
source('./R/plots.r')
source('./R/project.r')
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)

params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7,
                              sigma = 2, no_w=100)
object <- params 
sim <- project(params, t_max=1)

# extract n_pp and n from sim object 
n_pp <- sim@n_pp[1, ]
n <- sim@n[1, , ]
# sim@n[time,species, weight ]
# we need to get species index back even though there is only one species
dim(n) <- c(1, length(n))

mizerResult <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)

w0 <- object@w[1]

Beta <- log(object@species_params$beta)
sigma <- object@species_params$sigma
wFull <- params@w_full
xFull <- log(wFull)
dx <- xFull[2]-xFull[1]
s <- exp(-(xFull - xFull[1] - Beta)^2/(2*sigma^2))

#f <- object@dw_full*object@w_full*n_pp
f <- wFull*wFull*n_pp

energy <- dx*Re(fft(fft(s)*fft(f), inverse=TRUE)/length(s))
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
plot(log(object@w),energy[idx_sp])
lines(log(object@w), mizerResult)
# richard's code
params@interaction

summary(params)

params@pred_kernel

