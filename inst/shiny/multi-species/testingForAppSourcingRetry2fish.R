

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
l50 <-  c(1.0e+30,     2.0e+00,     5.0e+00,     6.0e+00,     3.0e+00,     3.6e+01,     8.0e+00,     5.1e+01) 
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
                                  h = 34, r_pp = 3*(10^6)*(10^(-2))))
# Choose kappa and weights properly
# distinguish between no_sp and no_extra_sp
# are the f0, h and lambda, r_pp chen here consistent with Mariella's choices ?
  
  # # #
  for (i in (1:no_sp)){
#for (i in (1:3)){
    a_m <- params_data$a2[i]
    b_m <- params_data$b2[i]
    L_inf_m <- params_data$Linf[i]
    L_mat <- params_data$Lmat[i]
    species_params <- data.frame(
      species = as.character(params_data$species[i]),
      w_min = params_data$Wegg[i],
      #w_inf = a_m*L_inf_m^b_m,
      w_inf = params_data$w_inf[i],
      #w_mat = a_m*L_mat^b_m,
      w_mat = params_data$w_mat[i],
      beta = params_data$beta[i],
      sigma = log(params_data$sigma[i]),
      z0 = 0,
      #alpha = input$alpha_anchovy, # unknown, mizer default=0.6
      #alpha = 0.6, # unknown, mizer default=0.6
      alpha = params@species_params$alpha[i],
      erepro = 0.1, # unknown, determined later
      sel_func = "sigmoid_length",
      gear = "sigmoid_gear",
      #l25 = l25["anchovy"],
      #l50 = l50["anchovy"],
      l25 = l25[i],
      l50 = l50[i],
      k = 0,
      k_vb = params_data$k_vb[i],
      a = a_m,
      b = b_m,
      #gamma = input$gamma_anchovy,
      #! does not look like we have a previous value for gamma or h
      #gamma = params_data$gamma[i],
      ##gamma = 0.00085,
      #!! when we turn off the gamma, the code currently goes wrong
      #! need to choose gamma better, 
      # gamma = params@species_params$gamma[i],
      #h = 50,
      ##h = params@species_params$h[i],
      linecolour = c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4", "darkorange", 
                     "darkorange1", "darkorange2", "darkorange3")[i],
      linetype = "solid"
    )
    #! SSB and effort need making dynamic
  #  p <- addSpecies(p, species_params, SSB = 2800, 
   #                 effort = 0, rfac = 1.01)
    #effort <- 1
    p <- addSpecies(p, species_params,
                    effort = effort, rfac = 1.01)
    
    
  }
  # # #
  
  sim <- project(p, t_max = 500, t_save = 5, effort = effort)
  #plot(sim)
  #########
  p <- sim@params
  no_sp <- length(p@species_params$species)
  no_t <- dim(sim@n)[1]
  p@initial_n <- sim@n[no_t, , ]
  p@initial_n_pp <- sim@n_pp[no_t, ]
  p@species_params$r_max <- Inf
  
  # Retune the values of erepro so that we get the correct level of
  # recruitment without stock-recruitment relationship
  mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = effort)
  gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
  rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
  for (i in (1:no_sp)) {
    gg0 <- gg[i, p@species_params$w_min_idx[i]]
    mumu0 <- mumu[i, p@species_params$w_min_idx[i]]
    DW <- p@dw[p@species_params$w_min_idx[i]]
    p@species_params$erepro[i] <- p@species_params$erepro[i] *
      (p@initial_n[i, p@species_params$w_min_idx[i]] *
         (gg0 + DW * mumu0)) / rdi[i]
  }
  
  
  sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
  plot(sim)
  
  # make effort variable


  p@species_params$erepro
  plotGrowthCurves(p, species=as.character(p@species_params$species[11:18]))
  plotGrowthCurves(p, species=as.character(p@species_params$species[12]))
  
  
  # Am running mariella's code and using csv to generate params information. Some values still need to be set, and 
  # we still need to implement fishing properly. 
  
  # decreased fishing effort to 0.6, and added growth curves
  
  
# Just pushed testingForAppSourcingRetry2fish.R This unfished system seems unstable, but is improved when fishing
  #is included
  # This code is like testingForAppSourcingRetry.R, but we have included proper l25 and l50 values, 
  # and are using mizer to determine the gamma. The system looks ok with fishing effort 0.6, but some erepro 
  # values are off, and there are some deviations in growth rates (although nothing really far out)
  
  # fixed l25 that was too big

  #made sure initialization runs for long enough to reach a steady state in fished system.