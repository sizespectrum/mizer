

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


####### recalculate pred_kernel directly from vignette form
beta1 <- object@species_params$beta
sigma1 <- object@species_params$sigma
pk1 = array(0,dim = c(dim(n)[1],dim(n)[2],length(n_pp)))
w1 <- object@w
wfull1 <- object@w_full

for (i in (1:dim(pk1)[1])){
  for (k in (1:dim(pk1)[2])){
    for (kp in (1:dim(pk1)[3])){
    
    if (wfull1[kp]<w1[k]) {
      pk1[i,k,kp] <- exp((-log(w1[k]/(wfull1[kp]*beta1[i]))^2)/(2*sigma1[i]))
    }
      }
    }
}




# Old get phi computer
# Inputs n, n_pp, object

n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * object@dw, "*", check.margin=FALSE) 
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
#phi_prey_species <- rowSums(sweep(object@pred_kernel[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
phi_prey_species <- rowSums(sweep(pk[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
#phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
phi_prey_background <- rowSums(sweep(pk,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
availableMaster <- phi_prey_species+phi_prey_background

### using new pk1

n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * object@dw, "*", check.margin=FALSE) 
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
#phi_prey_species <- rowSums(sweep(object@pred_kernel[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
phi_prey_species <- rowSums(sweep(pk1[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
#phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
phi_prey_background <- rowSums(sweep(pk1,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
availableMaster1 <- phi_prey_species+phi_prey_background

dim(availableMaster)

plot(availableMaster1[1,])

################### compare with fft based procedure ##########

availFFT <- getPhiPrey(object, n, n_pp)

# available energy using the old way, using the manually constructed pred kernel
plot(availableMaster1[1,])

# old code, old pred kernel
lines(availableMaster[1,])

# fft based new code for comparison
lines(availFFT[1,], col="Blue")


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

#### pred rate and getm2 calculator again for new pred_kernel
feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
n_total_in_size_bins <- sweep(n, 2, object@dw, '*', check.margin=FALSE) # N_i(w)dw
predd_rate <- sweep(pk1,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
m21 <- t(object@interaction) %*% colSums(aperm(predd_rate, c(2,1,3)),dims=1)[,idx_sp]


plot(m2[1,])
lines(m21[1,])


m2FFT <- getM2(object, n, n_pp)

lines(m2FFT[1,], col ="Blue")

# points, master using old pred_kernel
# solid blue, fft code
# solid black, master using new/manual pred_kernel

