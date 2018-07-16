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


all_efforts <- c(0, 1.4, 1.1)
names(all_efforts) <- c("knife_edge_gear", "sigmoid_gear", "sigmoid_gear_Anchovy")
effort <- all_efforts[1:2]
for (i in (1:no_sp)) {
    if (params_data$species[i] == "Anchovy") {
        effort <- c(effort, all_efforts[3])
    }
    a_m <- params_data$a2[i]
    b_m <- params_data$b2[i]
    L_inf_m <- params_data$Linf[i]
    L_mat <- params_data$Lmat[i]
    if (params_data$species[i] == "Anchovy") {
        gear <- "sigmoid_gear_Anchovy"
    } else {
        gear <- "sigmoid_gear"
    }
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
        gear = gear,
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
p <- steady(p, effort = effort, t_max = 100,  tol = 1e-3)

humboldt_params <- p
devtools::use_data(humboldt_params)

save(humboldt_params, file="humboldt_params.rda")
data("humboldt_params")
p <- humboldt_params

effort <- 1.4
sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)

p@species_params$erepro
plotGrowthCurves(p, species="Sardine")

# Now change one of the parameters

species_params <- p@species_params

species_params["Sardine", "gamma"] <- 
    species_params["Sardine", "gamma"] * 1.1
species_params["JMackerel", "w_mat"] <-
    species_params["JMackerel", "w_mat"] * 0.9
species_params["JMackerel", "w_min"] <-
    species_params["JMackerel", "w_min"] * 10

effort <- c(knife_edge_gear = 0, sigmoid_gear = 1.4, sigmoid_gear_Anchovy = 1.1)

pc <- MizerParams(
    species_params,
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
p <- steady(pc, effort = effort, t_max = 20, tol = 1e-2)

p@species_params$erepro
plotSpectra(p)


# Testing changing a general parameter

p_old <- p

p <- setBackground(
    set_scaling_model(
        #min_w_pp = input$min_w_pp, no_sp = input$no_bg_sp, no_w = 400,
        # min_w_pp = 1e-12, no_sp = input$no_bg_sp, no_w = 400,
        # min_w_inf = 2, max_w_inf = 6e5,
        # min_egg = 1e-4, min_w_mat = 2 / 10^0.6,
        # lambda = input$lambda, knife_edge_size = Inf,
        # f0 = input$f0, h = input$h_bkgd, r_pp = 10^input$log_r_pp
        min_w_pp = 1e-12,
        no_sp = 10, no_w = 400, min_w_inf = 2, max_w_inf = 6e5,
        min_egg = 1e-4, min_w_mat = 2 / 10^0.6, 
        lambda = 2.12,
        knife_edge_size = Inf,
        
    )
)
# Loop over all foreground species and add them one-by-one to the new
# background
effort <- 0
names(effort) <- "knife_edge_gear"
all_efforts <- c(0, 1.4, 1.1)
names(all_efforts) <- c("knife_edge_gear", "sigmoid_gear", "sigmoid_gear_Anchovy")
no_sp <- length(p_old@A)
for (sp in (1:no_sp)[!is.na(p_old@A)]) {
    if (!(p_old@species_params[sp, "gear"] %in% names(effort))) {
        effort <- c(effort, all_efforts[p_old@species_params[sp, "gear"]])
    }
    p <- addSpecies(p, p_old@species_params[sp, ],
                    effort = effort,
                    rfac=Inf)
}

# Run to steady state
p <- steady(p, effort = effort, 
            t_max = 100, tol = 1e-2)

############################################################################

plotGrowthCurves(p, species="Sardine")


getAgesAtMaturity <- function(p){
  no_sp <- length(p@species_params$species)
  age_at_maturity <- rep(2,no_sp)
  for (i in (1:no_sp)[!is.na(p@A)]){
    age_at_maturity[i] <- -log(1-((p@species_params$w_mat[i]/p@species_params$w_inf[i])^(1/p@species_params$b[i]))
)/p@species_params$k_vb[i]  
  }
 # age_at_maturity[is.na(p@A)] <- 2
  return(age_at_maturity)
}


age_at_maturity <- getAgesAtMaturity(p)

#((p@species_params$w_mat[i]/p@species_params$w_inf[i])^(1/p@species_params$b[i]))

# actually use VB to compute this
retune_h_for_growth <- function(p,age_at_maturity){
  q <- p
  for (i in  (1:no_sp)[!is.na(p@A)]){
    q@species_params$h[i] <- p@species_params$ks[i] +
      (p@species_params$w_mat[i]^(1-p@n) - p@species_params$w_min[i]^(1-p@n))/(
        age_at_maturity[i]*(1-p@n)*p@species_params$alpha[i]*p@f0
      )
  }
  return(q)
}


qq <- retune_h_for_growth(p,age_at_maturity)
plotGrowthCurves(qq, species="Sardine")
qq <- steady(qq, effort = effort, 
            t_max = 100, tol = 1e-2)
plotGrowthCurves(qq, species="Sardine")


retune_h_and_ks_for_growth <- function(p,age_at_maturity){
  q <- p
  for (i in  (1:no_sp)[!is.na(p@A)]){
    q@species_params$h[i] <- (p@species_params$w_mat[i]^(1-p@n) - p@species_params$w_min[i]^(1-p@n))/(
        age_at_maturity[i]*(1-p@n)*p@species_params$alpha[i]*p@f0 -0.2
      )
    q@species_params$k_s[i] <-0.2* q@species_params$h[i] 
  }
  return(q)
}

