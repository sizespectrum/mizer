


####################################### run old code up to break point




source("C:/Users/user/Desktop/mizer projects/retune/mizer/inst/shiny/multi-species/ModeloMIZER_NCME_v1DirAlt.R")


effort <- 0.6

params_data <- read.csv("C:/Users/user/Desktop/mizer projects/retune/mizer/inst/shiny/multi-species/speciesNCME_edited2.csv")

no_sp <- dim(params_data)[1]
#l50 <- c(12.5, 15.48, 15.48)
l50 <- rep(12.5,no_sp)
#names(l50) <- c("anchovy", "jmackerel", "jmackerelB")
names(l50) <- as.character(params_data$species)
#sd <- c(0.462, 2.10, 2.10)
sd <- rep(2.10,no_sp)
l25 = l50 - log(3) * sd


# changed l25[5] from 3 to 2.9
# changed l50[2] from 2 to 1.9
# changed l25[7] to from 4.9e+01 to 4.9e+00
l25 <- c(1.0e+29,     1.9e+00,     4.0e+00,     5.0e+00,     2.9e+00,     3.2e+01,     4.9e+00,     4.9e+01 )
l50 <-  c(1.1e+29,     2.0e+00,     5.0e+00,     6.0e+00,     3.0e+00,     3.6e+01,     8.0e+00,     5.1e+01) 
names(l25) <- as.character(params_data$species)
names(l50) <- as.character(params_data$species)


# do we want to make ks variable too ?

p <- setBackground(set_scaling_model(no_sp = 10, no_w = 400,
                                     min_w_inf = 2, max_w_inf = 6e5,
                                     min_egg = 1e-4, min_w_mat = (1/5)* 10^(0.4),
                                     #knife_edge_size = Inf, kappa = 10000,
                                     knife_edge_size = Inf, kappa = params@kappa,
                                     #lambda = 2.08, f0 = 0.6,
                                     lambda = 1.95, f0 = params@f0,
                                     h = 34, r_pp = (10^(-2))))
#########################


#i <- 1
for (ij in (1:no_sp)){
  
#for (i in (1:3)){
a_m <- params_data$a2[ij]
b_m <- params_data$b2[ij]
L_inf_m <- params_data$Linf[ij]
L_mat <- params_data$Lmat[ij]
species_params <- data.frame(
  species = as.character(params_data$species[ij]),
  w_min = params_data$Wegg[ij],
  #w_inf = a_m*L_inf_m^b_m,
  w_inf = params_data$w_inf[ij],
  #w_mat = a_m*L_mat^b_m,
  w_mat = params_data$w_mat[ij],
  beta = params_data$beta[ij],
  sigma = log(params_data$sigma[ij]),
  z0 = 0,
  #alpha = input$alpha_anchovy, # unknown, mizer default=0.6
  #alpha = 0.6, # unknown, mizer default=0.6
  alpha = params@species_params$alpha[ij],
  erepro = 0.1, # unknown, determined later
  sel_func = "sigmoid_length",
  gear = "sigmoid_gear",
  #l25 = l25["anchovy"],
  #l50 = l50["anchovy"],
  l25 = l25[ij],
  l50 = l50[ij],
  k = 0,
  k_vb = params_data$k_vb[ij],
  a = a_m,
  b = b_m
  #gamma = input$gamma_anchovy,
  #! does not look like we have a previous value for gamma or h
  #gamma = params_data$gamma[i],
  ##gamma = 0.00085,
  #!! when we turn off the gamma, the code currently goes wrong
  #! need to choose gamma better, 
  # gamma = params@species_params$gamma[i],
  #h = 50,
  ##h = params@species_params$h[i],
  #linecolour = c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4", "darkorange", 
  #            "darkorange1", "darkorange2", "darkorange3")[i],
  #linetype = "solid"
)
#! SSB and effort need making dynamic
#  p <- addSpecies(p, species_params, SSB = 2800, 
#                 effort = 0, rfac = 1.01)
#effort <- 1

################################################ define retune abundance



retune_abundance <- function(params, retune) {
  no_sp <- length(params@species_params$w_inf)  # Number of species
  if (length(retune) != no_sp) {
    stop("retune argument has the wrong length")
  }
  # We try to match the original abundance between the maturity size
  # of the smallest species and the maximum size of the largest species.
  # Determine the indices of these limits
  idx_start <- sum(params@w <= min(params@species_params$w_mat))
  idx_stop <- sum(params@w < max(params@species_params$w_inf))
  # More precisely, we find the abundance multipliers A_i so
  # that the integral of the square of the relative distance 
  # (sum_{i not in L} A_i*N_i(w) + sum_{i not in L} N_i(w) - sc(w))/sc(w) 
  # over w, between our limits, is minimized, where  L is the set of all
  # retuneable species.
  
  # deal with zero entries in params@sc
  nonzero <- params@sc > 0
  sc <- params@sc[nonzero]
  # rho is the total abundance of all the non-tunable species
  rho <- colSums(params@initial_n[!retune, nonzero, drop=FALSE])
  
  # Use Singular Value Decomposition to find optimal abundance multipliers.
  # See Numerical Recipes section 15.4.2
  #
  # Rescale by sc
  A <- t(sweep(params@initial_n[retune, nonzero, drop=FALSE], 2, sc, "/"))
  b <- (sc - rho) / sc
  
  sv <- svd(A)
  di <- 1/sv$d  # inverse of singular values
  di[di > 10^8] <- 0  # cut off
  x <- sweep(sv$v, 2, di, "*") %*% t(sv$u) %*% b
  A2 <- rep(1, no_sp) 
  A2[retune] <- x
  # We may have to repeat this if any of the multipliers is negative
  if (any(A2 <= 0)) {
    # Set abundance of those species to tiny
    params@initial_n[A2 <= 0, ] <- max(params@initial_n[A2 <= 0, ] * 1e-40, 1e-40)
    # and try again retuning the remaining retunable species
    retune <- retune & (A2 > 0)
    if (any(retune)) {
      params <- retune_abundance(params, retune)
    } else {
      stop("Unable to retune.")
    }
  } else {
    # Use these abundance multipliers to rescale the abundance curves
    params@initial_n <- params@initial_n * A2
    # update SSB
    params@A <- params@A * A2
  }
  return(params)
}

################################ get inputs for the following add species internal code, and test

#p <- addSpecies(p, species_params,
#                effort = effort, rfac = 1.01)
params <- p
rfac <- 1.01
#SSB <- NA
SSB <- 85900583
##################################

if (any(species_params$species %in% params@species_params$species)) {
  stop("You can not add species that are already there.")
}
# calculate h if it is missing
if (!hasName(species_params, "h") || is.na(species_params$h)) {
  message("Note: \tNo h column in new species data frame so using f0 and k_vb to
          calculate it.")
  if(!hasName(species_params, "k_vb")) {
    stop("\t\tExcept I can't because there is no k_vb column in the new species data frame")
  }
  fc <- 0.2/species_params$alpha
  species_params$h <- 3*species_params$k_vb*(species_params$w_inf^(1/3))/
    (species_params$alpha*params@f0-0.2)
}

# calculate ks if it is missing
if (!hasName(species_params, "ks") || is.na(species_params$ks)){
  message("Note: \tNo ks column in new species data frame. Setting ks = 0.2*h.")
  species_params$ks <- 0.2*species_params$h # mizer's default setting
}

# calculate gamma if it is missing
if (!hasName(species_params, "gamma") || is.na(species_params$gamma)){
  message("Note: \tNo gamma column in new species data frame so using f0, h, beta, sigma, lambda and kappa to calculate it.")
  ae <- sqrt(2*pi) * species_params$sigma * 
    species_params$beta^(params@lambda-2) * 
    exp((params@lambda-2)^2 * species_params$sigma^2 / 2)
  species_params$gamma <- (species_params$h / (params@kappa * ae)) * 
    (params@f0 / (1 - params@f0))
}

# calculate w_min_idx
species_params$w_min_idx <- sum(params@w<=species_params$w_min)

# provide erepro column that is later overwritten
species_params$erepro <- 0.1

# Move linecolour and linetype into species_params
params@species_params$linetype <- 
  params@linetype[params@species_params$species]
params@species_params$linecolour <- 
  params@linecolour[params@species_params$species]
# TODO: Check if we need to do this with selectivity as well

# Make sure that all columns exist in both data frames
missing <- setdiff(names(params@species_params), names(species_params))
species_params[missing] <- NA
missing <- setdiff(names(species_params), names(params@species_params))
params@species_params[missing] <- NA

# add the new species (with parameters described by species_params), 
# to make a larger species_params dataframe.
combi_species_params <- rbind(params@species_params, species_params)

# use dataframe and global settings from params to make a new MizerParams 
# object.
p <- MizerParams(
  combi_species_params,
  p = params@p,
  n = params@n,
  q = params@q,
  lambda = params@lambda,
  f0 = params@f0,
  kappa = params@kappa,
  min_w = min(params@w),
  max_w = max(params@w),
  no_w = length(params@w),
  min_w_pp = min(params@w_full),
  w_pp_cutoff = max(params@w_full),
  r_pp = (params@rr_pp / (params@w_full ^ (params@p - 1)))[1]
)
# Use the same resource spectrum as params
p@initial_n_pp <- params@initial_n_pp
p@cc_pp <- params@cc_pp
new_sp <- length(params@species_params$species) + 1
no_sp <- new_sp
# Initially use abundance curves for pre-existing species 
# (we shall retune the abundance multipliers of such 
# species from the background later)
p@initial_n[1:(new_sp - 1), ] <- params@initial_n
# Use the same psi and mu_b as before for old species
p@psi[1:(new_sp - 1), ] <- params@psi
p@sc <- params@sc
p@mu_b[1:(new_sp - 1), ] <- params@mu_b
p@mu_b[new_sp, ] <- params@mu_b[1, ]
p@srr <- params@srr

# Turn off self-interaction of the new species, so we can determine the
# growth rates, and death rates induced upon it by the pre-existing species
p@interaction[new_sp, new_sp] <- 0
# compute death rate for new species
mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = effort)[new_sp, ]
# compute growth rate for new species
gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)[new_sp, ]

# Compute solution for new species
w_inf_idx <- sum(p@w < p@species_params$w_inf[new_sp])
idx <- p@species_params$w_min_idx[new_sp]:(w_inf_idx-1)
if (any(gg[idx]==0)) {
  stop("Can not compute steady state due to zero growth rates")
}
p@initial_n[new_sp, ] <- 0
p@initial_n[new_sp, p@species_params$w_min_idx[new_sp]:w_inf_idx] <- 
  c(1, cumprod(gg[idx] / ((gg + mumu * p@dw)[idx+1])))
if (any(is.infinite(p@initial_n))) {
  stop("Candidate steady state holds infinities")
}
if (any(is.na(p@initial_n) || is.nan(p@initial_n))) {
  stop("Candidate steady state holds none numeric values")
}

# Normalise solution
if (is.na(SSB)) {
  # If spawning stock biomass of new species is not supplied, 
  # normalise solution so that at its maximum it lies at half the 
  # power law, and then calculate its SSB.
  # We choose the maximum of the biomass density in log space
  # because that is always an increasing function at small size.
  idx <- which.max(p@initial_n[new_sp, ] * p@w^p@lambda)
  p@initial_n[new_sp, ] <- p@initial_n[new_sp, ] *
    p@kappa * p@w[idx]^(-p@lambda) / p@initial_n[new_sp, idx] / 2
  SSB <- sum(p@initial_n[new_sp, ] * p@w * p@dw * p@psi[new_sp, ])
} else {
  unnormalised_SSB <- sum(p@initial_n[new_sp,] * p@w * p@dw * 
                            p@psi[new_sp, ])
  p@initial_n[new_sp, ] <- p@initial_n[new_sp, ] * SSB / unnormalised_SSB
}
p@A <- c(params@A, SSB)

# Turn self interaction back on
p@interaction[new_sp, new_sp] <- 1

# Retune the abundance multipliers to recreate the aggregate abundance
# spectrum of the old params object.
# First identify the retunable species. These are all background
# species except the largest one
retune <- is.na(p@A)
largest_back_idx <- which.max(p@species_params$w_inf[retune])
retune[largest_back_idx] <- FALSE
p <- retune_abundance(p, retune)

# Retune the values of erepro, so that we are at steady state.
# First get death and growth rates
mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = effort)
gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
erepro_final <- 1:no_sp  # set up vector of right dimension
for (i in (1:no_sp)) {
  gg0 <- gg[i, p@species_params$w_min_idx[i]]
  mumu0 <- mumu[i, p@species_params$w_min_idx[i]]
  DW <- p@dw[p@species_params$w_min_idx[i]]
  if (!rdi[i]==0){
    erepro_final[i] <- p@species_params$erepro[i] *
      (p@initial_n[i, p@species_params$w_min_idx[i]] *
         (gg0 + DW * mumu0)) / rdi[i]
  }
  else {
    erepro_final[i] <- 0.1
  }
}
if (is.finite(rfac)) {
  # erepro has been multiplied by a factor of (rfac/(rfac-1)) to
  # compensate for using a stock recruitment relationship.
  erepro_final <- (rfac / (rfac - 1)) * erepro_final
}
p@species_params$erepro <- erepro_final

p@species_params$r_max <- p@species_params$w_inf
# set rmax=fac*RDD
# note that erepro has been multiplied by a factor of (rfac/(rfac-1)) to
# compensate for using a stock recruitment relationship.
p@species_params$r_max <-
  (rfac - 1) * getRDI(p, p@initial_n, p@initial_n_pp)[,1]
#return(p)
#}

sim <- project(p, t_max = 50, t_save = 5, effort = effort, dt=0.1)
#plot(sim)
#########
p <- sim@params
no_sp <- length(p@species_params$species)
no_t <- dim(sim@n)[1]
p@initial_n <- sim@n[no_t, , ]
p@initial_n_pp <- sim@n_pp[no_t, ]
#p@species_params$r_max <- Inf

# Retune foreground species to have the correct abundances
foreground <- (1:no_sp)[!is.na(p@A)]
for (sp in foreground) {
  SSB <- p@A[sp]
  unnormalised_SSB <- sum(p@initial_n[sp,] * p@w * p@dw * 
                            p@psi[sp, ])
  p@initial_n[sp, ] <- p@initial_n[sp, ] * SSB / unnormalised_SSB
}

# Retune the values of erepro so that we get the correct level of
# recruitment without stock-recruitment relationship
mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = effort)
gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
for (ii in (1:no_sp)) {
  gg0 <- gg[ii, p@species_params$w_min_idx[ii]]
  mumu0 <- mumu[i, p@species_params$w_min_idx[ii]]
  DW <- p@dw[p@species_params$w_min_idx[ii]]
  p@species_params$erepro[i] <- p@species_params$erepro[ii] *
    (p@initial_n[i, p@species_params$w_min_idx[ii]] *
       (gg0 + DW * mumu0)) / rdi[ii]
}

}
p@species_params$r_max <- Inf

sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)

