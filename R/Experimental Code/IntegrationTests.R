

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


####### recalculate pred_kernel directly from vignette form
beta1 <- object@species_params$beta
sigma1 <- object@species_params$sigma
pk1 = array(0,dim = c(dim(n)[1],dim(n)[2],length(n_pp)))
w1 <- object@w
wfull1 <- object@w_full

for (i in (1:dim(pk1)[1])){
  for (k in (1:dim(pk1)[2])){
    for (kp in (1:dim(pk1)[3])){
    x1 <- log(w1)
    xfull1 <- log(wfull1)
    
    #if (wfull1[kp]<=w1[k]) {
    if (-log(beta1[i])-3*sigma1[i]<= xfull1[kp]-x1[k]) {
      if (xfull1[kp]-x1[k]<=0) {
      
    
      pk1[i,k,kp] <- exp((-log(w1[k]/(wfull1[kp]*beta1[i]))^2)/(2*sigma1[i]^2))
    }}
      }
    }
}

# Old get phi computer (with new widths)

xfull<- log(object@w_full)
newWidth<-object@w_full*(xfull[2]-xfull[1])
newWidthshort<-object@w*(xfull[2]-xfull[1])
#n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * object@dw, "*", check.margin=FALSE) 
n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * newWidthshort, "*", check.margin=FALSE) 
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
#phi_prey_species <- rowSums(sweep(object@pred_kernel[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
phi_prey_species <- rowSums(sweep(pk1[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
#phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
#phi_prey_background <- rowSums(sweep(pk1,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
phi_prey_background <- rowSums(sweep(pk1,3,object@w_full*newWidth*n_pp,"*", check.margin=FALSE),dims=2)
availableMasterNW <- phi_prey_species+phi_prey_background

# Old get phi computer (with old widths)

xfull<- log(object@w_full)
newWidth<-object@w_full*(xfull[2]-xfull[1])
newWidthshort<-object@w*(xfull[2]-xfull[1])
n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * object@dw, "*", check.margin=FALSE) 
#n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * newWidthshort, "*", check.margin=FALSE) 
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
#phi_prey_species <- rowSums(sweep(object@pred_kernel[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
phi_prey_species <- rowSums(sweep(pk1[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
#phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
phi_prey_background <- rowSums(sweep(pk1,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
phi_prey_background <- rowSums(sweep(pk1,3,object@w_full*newWidth*n_pp,"*", check.margin=FALSE),dims=2)
availableMasterOW <- phi_prey_species+phi_prey_background

availFFT <- getPhiPrey(object, n, n_pp)


plot(availableMasterOW[1,])
lines(availableMasterNW[1,], col="Red")
# fft based new code for comparison
lines(availFFT[1,], col="Blue")

##############################

##### old pred rate calculator

feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
n_total_in_size_bins <- sweep(n, 2, object@dw, '*', check.margin=FALSE) # N_i(w)dw
#predd_rate <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)
predd_rate <- sweep(pk,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)

#### old getM2 calculator
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
m2OW <- t(object@interaction) %*% colSums(aperm(predd_rate, c(2,1,3)),dims=1)[,idx_sp]

xfull<- log(object@w_full)
newWidth<-object@w_full*(xfull[2]-xfull[1])
newWidthshort<-object@w*(xfull[2]-xfull[1])


#### pred rate and getm2 calculator again for new pred_kernel
feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp)
#n_total_in_size_bins <- sweep(n, 2, object@dw, '*', check.margin=FALSE) # N_i(w)dw
n_total_in_size_bins <- sweep(n, 2, newWidthshort, '*', check.margin=FALSE) # N_i(w)dw
predd_rate <- sweep(pk1,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
m2NW <- t(object@interaction) %*% colSums(aperm(predd_rate, c(2,1,3)),dims=1)[,idx_sp]

m2FFT <- getM2(object, n, n_pp)


plot(m2OW[1,])
lines(m2NW[1,], col="Red")

lines(m2FFT[1,], col ="Blue")

# points, master using old pred_kernel
# solid red, master with new spacing
# blue fft result

