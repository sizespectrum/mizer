#####################
# Source the mizer code without loading it as a package
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

params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7)
object <- params 
sim <- project(params, t_max=1)

# extract n_pp and n from sim object 
n_pp <- sim@n_pp[1, ]
n <- sim@n[1, , ]

# sim@n[time,species, weight ]
# we need to get species index back even though there is only one species
dim(n) <- c(1, length(n))

params@interaction

summary(params)

params@pred_kernel

#####################

# This is the value you want to reproduce with fft method:
phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)

object@species_params$beta
object@species_params$sigma
object@w[1]
log(object@dw[1])

log(object@w[2])-log(object@w[1])


wFull <- params@w_full 


w0 <- eggsize

# how to find eggsize ?

xFull <- log(wFull/w0)

# how to find sigma and Beta ? note, our Beta = log(beta_from_vignette) ?

s <- exp(-(xFull -Beta)^2/(2*sigma^2))

N <- length(s)

dx <- xFull[2]-xFull[1]

# n[j,wp]= N_j(w_p)

f <- sweep(sweep((object@interaction %*% n), 1, n_pp, "+"), 2, wFull, "*")[i]

# could it be f <- sweep(sweep((object@interaction %*% n), 1, n_pp, "+"), 2, wFull^2, "*")[i] ?

# Alternative loop based approach...

A <- (object@interaction %*% n)
for (i in 1:(dim(object@interaction)[1])){
 A[i,] <- A[i,] + nn_p}

for (i in 1:(dim(object@interaction)[2])){
    A[,i] <- A[,i] * wFull}

f <- A[i, ]


# can/should we do the convolution integral for all species i simultaineously ?

energy <- dx*Re(fft(fft(s)*fft(f), inverse=TRUE)/N)

