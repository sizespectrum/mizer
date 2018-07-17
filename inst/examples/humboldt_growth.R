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

plotGrowthCurves(p, species="Sardine")


#################################################### function to return age[i,w]
getMizerGrowthCurve <- function(p, max_age = 40){
  age <- seq(0, max_age, length.out = 50)
  ages_at_weight <- p@psi
  ages_at_weight[] <- 0
  no_sp <- length(p@species_params$species)
  for (i in 1:no_sp){
    g <- getEGrowth(p, p@initial_n, p@initial_n_pp)
    
    
    g_fn <- stats::approxfun(p@w, g[i, ])
    myodefun <- function(t, state, parameters){
      return(list(g_fn(state)))
    }
    weights <- deSolve::ode(y = p@species_params$w_min[i], 
                            times = age, func = myodefun)[,2]
    weights_ref <- weights
    #weights_ref[p@w>p@species_params$w_inf[i]] <- 0
    fillinages <- stats::approxfun(weights_ref, age)
    ages_at_weight[i,] <- fillinages(p@w)
    
  }
  return(ages_at_weight)}

p@species_params$species[9]

################ compare function to growth curve of sardine
HH <- getMizerGrowthCurve(p)[9,]
  
  HH[is.na(HH)] <- 0

plot(HH,p@w,log="y")
abline(h=p@species_params$w_mat[9])
plotGrowthCurves(p, species="Sardine")

############## determine ages at maturity in mizer (note slight discrepencies due to 
# some species characteristic weights not lying on grid points )

ages_at_maturity_mizer <- function(p, max_age = 400){
  data <- getMizerGrowthCurve(p, max_age)
  no_sp <- length(p@species_params$species)
  output <- p@species_params$w_mat
  for (i in 1:no_sp){
  output[i] <- data[i,sum(p@w<=p@species_params$w_mat[i])]
  }
return(output)
}
ages_at_maturity_mizer(p)


############### give approximate theoretical approximations to ages at maturity

ages_at_maturity_theory <- function(p){
  return(
    (p@species_params$w_mat^(1-p@n) - p@species_params$w_min^(1-p@n))/(
      (1-p@n)*(p@species_params$alpha*p@f0*p@species_params$h - p@species_params$ks) 
    )
  )
  
}
ages_at_maturity_theory(p)

# there is some deviation, so I compute ages in a second way to avoid issues about grid points

ages_at_maturity_mizer2 <- function(p, max_age = 400){
  age <- seq(0, max_age, length.out = 50)
  ages_at_weight <- p@species_params$w_mat
  no_sp <- length(p@species_params$species)
  for (i in 1:no_sp){
    g <- getEGrowth(p, p@initial_n, p@initial_n_pp)
    
    
    g_fn <- stats::approxfun(p@w, g[i, ])
    myodefun <- function(t, state, parameters){
      return(list(g_fn(state)))
    }
    weights <- deSolve::ode(y = p@species_params$w_min[i], 
                            times = age, func = myodefun)[,2]
    weights_ref <- weights
    #weights_ref[p@w>p@species_params$w_inf[i]] <- 0
    fillinages <- stats::approxfun(weights_ref, age)
    ages_at_weight[i] <- fillinages(p@species_params$w_mat[i])
    
  }
  return(ages_at_weight)}

ages_at_maturity_mizer2(p)
ages_at_maturity_theory(p)
# there is quite some deviation between these systems.


############################ compute VB growth curve, 

################################# invert growth, and get target ages at maturuty

################################# demonstrait that we can retune h to fit

############################### make new code with ks and h inititalization
