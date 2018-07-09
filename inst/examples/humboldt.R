params_data <- read.csv(system.file("extdata", "speciesNCME_edited2.csv", package = "mizer"))
effort <- 0.6

no_sp <- dim(params_data)[1]

l25 <- c(1.0e+29,     1.9e+00,     4.0e+00,     5.0e+00,     2.9e+00,     3.2e+01,     4.9e+00,     4.9e+01 )
l50 <-  c(1.1e+29,     2.0e+00,     5.0e+00,     6.0e+00,     3.0e+00,     3.6e+01,     8.0e+00,     5.1e+01) 
names(l25) <- as.character(params_data$species)
names(l50) <- as.character(params_data$species)

p <- setBackground(
    set_scaling_model(
        no_sp = 10, no_w = 400, min_w_inf = 2, max_w_inf = 6e5,
        min_egg = 1e-4, min_w_mat = 2 / 10^0.6, 
        knife_edge_size = Inf)
)
plotSpectra(p)

pa <- c()
for (i in (1:no_sp)) {
    a_m <- params_data$a2[i]
    b_m <- params_data$b2[i]
    L_inf_m <- params_data$Linf[i]
    L_mat <- params_data$Lmat[i]
    species_params <- data.frame(
        species = as.character(params_data$species[i]),
        w_min = params_data$Wegg[i],
        w_inf = params_data$w_inf[i],
        w_mat = params_data$w_mat[i],
        beta = params_data$beta[i],
        sigma = log(params_data$sigma[i]),
        z0 = 0,
        alpha = 0.6,
        erepro = 0.1, # unknown, determined later
        sel_func = "sigmoid_length",
        gear = "sigmoid_gear",
        l25 = l25[i],
        l50 = l50[i],
        k = 0,
        k_vb = params_data$k_vb[i],
        a = a_m,
        b = b_m
    )

    p <- addSpecies(p, species_params, effort = effort)
    
    # Run to steady state
    t_max = 50
    # Force the recruitment to stay at this level
    rdd <- getRDD(p, p@initial_n, p@initial_n_pp)
    p@srr <- function(rdi, species_params) {rdd}
    sim_old <- project(p, t_max = t_max, t_save = 50, effort = effort)
    # Restore original stock-recruitment relationship
    p@srr <- params@srr
    plotBiomass(sim_old)
    
    no_sp <- length(p@species_params$species)
    no_t <- dim(sim_old@n)[1]
    p@initial_n <- sim_old@n[no_t, , ]
    p@initial_n_pp <- sim_old@n_pp[no_t, ]
    
    # Retune the values of erepro so that we get the correct level of
    # recruitment
    mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = effort)
    gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
    rdd <- getRDD(p, p@initial_n, p@initial_n_pp)
    # TODO: vectorise this
    for (i in (1:no_sp)) {
        gg0 <- gg[i, p@species_params$w_min_idx[i]]
        mumu0 <- mumu[i, p@species_params$w_min_idx[i]]
        DW <- p@dw[p@species_params$w_min_idx[i]]
        p@species_params$erepro[i] <- p@species_params$erepro[i] *
            (p@initial_n[i, p@species_params$w_min_idx[i]] *
                 (gg0 + DW * mumu0)) / rdd[i]
    }
    #p <- steady(p, effort = effort, t_max = 10)
    pa <- c(p, pa)
}

sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)

p@species_params$erepro
plotGrowthCurves(p, species=as.character(p@species_params$species[11:18]))
plotGrowthCurves(p, species=as.character(p@species_params$species[12]))
  