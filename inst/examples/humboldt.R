params_data <- read.csv(system.file("extdata", "speciesNCME_edited2.csv", package = "mizer"))
effort <- 1.4

no_sp <- dim(params_data)[1]

l25 <- c(1.0e+29,     1.9e+00,     4.0e+00,     5.0e+00,     2.9e+00,     3.2e+01,     4.9e+00,     4.9e+01 )
l50 <-  c(1.1e+29,     2.0e+00,     5.0e+00,     6.0e+00,     3.0e+00,     3.6e+01,     8.0e+00,     5.1e+01) 
names(l25) <- as.character(params_data$species)
names(l50) <- as.character(params_data$species)

p <- setBackground(
    set_scaling_model(min_w_pp = 1e-12,
        no_sp = 10, no_w = 400, min_w_inf = 2, max_w_inf = 6e5,
        min_egg = 1e-4, min_w_mat = 2 / 10^0.6, 
        lambda = 2.12,
        knife_edge_size = Inf)
)

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

    p <- addSpecies(p, species_params, effort = effort, rfac=Inf)
}

# Run to steady state
p <- steady(p, effort = effort, t_max = 100,  tol = 1e-2)

sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)

p@species_params$erepro
plotGrowthCurves(p, species="Sardine")

# Now change one of the parameters

p@species_params["Sardine", "gamma"] <- 
    p@species_params["Sardine", "gamma"] * 1.1

pc <- MizerParams(
    p@species_params,
    p = p@p,
    n = p@n,
    q = p@q,
    lambda = p@lambda,
    f0 = p@f0,
    kappa = p@kappa,
    min_w = min(p@w),
    max_w = max(p@w),
    no_w = length(p@w),
    min_w_pp = min(p@w_full),
    w_pp_cutoff = max(p@w_full),
    r_pp = (p@rr_pp / (p@w_full ^ (p@p - 1)))[1]
)
pc@linetype <- p@linetype
pc@linecolour <- p@linecolour
pc@A <- p@A
pc@sc <- p@sc
pc@cc_pp <- p@cc_pp
pc@mu_b <- p@mu_b

pc@initial_n <- p@initial_n
pc@initial_n_pp <- p@initial_n_pp
# Run to steady state
#p <- steady(pc, effort = effort, t_max = 10, plot = TRUE)
p <- steady(pc, effort = effort, t_max = 10)
