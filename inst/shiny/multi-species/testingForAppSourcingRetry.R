
source("C:/Users/user/Desktop/mizer projects/retune/mizer/inst/shiny/multi-species/ModeloMIZER_NCME_v1DirAlt.R")

params_data <- read.csv("C:/Users/user/Desktop/mizer projects/scaling/mizer/inst/shiny/multi-species/speciesNCME_edited2.csv")

no_sp <- dim(params_data)[1]
#l50 <- c(12.5, 15.48, 15.48)
l50 <- rep(12.5,no_sp)
#names(l50) <- c("anchovy", "jmackerel", "jmackerelB")
names(l50) <- as.character(params_data$species)
#sd <- c(0.462, 2.10, 2.10)
sd <- rep(2.10,no_sp)
l25 = l50 - log(3) * sd


# do we want to make ks variable too ?

p <- setBackground(set_scaling_model(no_sp = 10, no_w = 400,
                                  min_w_inf = 10, max_w_inf = 1e5,
                                  min_egg = 1e-4, min_w_mat = 10^(0.4),
                                  knife_edge_size = Inf, kappa = 10000,
                                  #lambda = 2.08, f0 = 0.6,
                                  lambda = params@lambda, f0 = params@f0,
                                  h = 34, r_pp = 10^(-2)))
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
      sigma = params_data$sigma[i],
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
      gamma = 0.00085,
      #! need to choose gamma better, 
      # gamma = params@species_params$gamma[i],
      #h = 50,
      h = params@species_params$h[i],
      linecolour = c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4", "darkorange", 
                     "darkorange1", "darkorange2", "darkorange3")[i],
      linetype = "solid"
    )
    #! SSB and effort need making dynamic
    p <- addSpecies(p, species_params, SSB = 2800, 
                    effort = 0, rfac = 1.01)
    
  }
  # # #
  
  sim <- project(p, t_max = 50, t_save = 5, effort = 0)
  
  #########
  p <- sim@params
  no_sp <- length(p@species_params$species)
  no_t <- dim(sim@n)[1]
  p@initial_n <- sim@n[no_t, , ]
  p@initial_n_pp <- sim@n_pp[no_t, ]
  p@species_params$r_max <- Inf
  
  # Retune the values of erepro so that we get the correct level of
  # recruitment without stock-recruitment relationship
  mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = 0)
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
  
  
  sim <- project(p, t_max = 15, t_save = 0.1, effort = 0)
  plot(sim)
  
  # make effort variable


# Am running mariella's code and using csv to generate params information. Some values still need to be set, and 
  # we still need to implement fishing properly. 
