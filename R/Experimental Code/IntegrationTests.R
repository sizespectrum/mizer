

library(mizer)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)

data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
object <- params
n <- get_initial_n(object)
n_pp <- object@cc_pp

######## calculate pred_kernel in old style ##########
pk = array(0,dim = c(dim(n)[1],dim(n)[2],length(n_pp)))
pk[] <- object@species_params$beta
pk <- exp(-0.5*sweep(log(sweep(sweep(pk,3,object@w_full,"*")^-1,2,object@w,"*")),1,object@species_params$sigma,"/")^2)
combn <- NULL
pk <- sweep(pk,c(2,3),combn(object@w_full,1,function(x,w)x<w,w=object@w),"*") # find out the untrues and then multiply


#######


# Old get phi computer
# Inputs n, n_pp, object

n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * object@dw, "*", check.margin=FALSE) 
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
#phi_prey_species <- rowSums(sweep(object@pred_kernel[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
phi_prey_species <- rowSums(sweep(pk[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
#phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
phi_prey_background <- rowSums(sweep(pk,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
availableMaster <- phi_prey_species+phi_prey_background

dim(availableMaster)

plot(availableMaster[1,])

################### compare with fft based procedure ##########

availFFT <- getPhiPrey(object, n, n_pp)
lines(availFFT[1,])


# can compare this with fft available energy, also test if the absence of truncated 
# feeding kernel makes any difference.

##### old pred rate calculator

feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
n_total_in_size_bins <- sweep(n, 2, object@dw, '*', check.margin=FALSE) # N_i(w)dw
#predd_rate <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)
predd_rate <- sweep(pk,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)

#### old getM2 calculator
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
m2 <- t(object@interaction) %*% colSums(aperm(predd_rate, c(2,1,3)),dims=1)[,idx_sp]

dim(m2)
plot(m2[1,])

m2FFT <- getM2(object, n, n_pp)

lines(m2FFT[1,])

