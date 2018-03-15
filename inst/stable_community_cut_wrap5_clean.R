
#19 prepared (in stable_community_cut_wrap5_clean.R) cleaner version 
#of set_scaling_model, and 
# pasted it into wrapper_functions.R. It seems like the resulting 
# mizer build ok, but for some reason set_scaling_model() does not 
# work after that. However, if one runs wrapper_functions.R 
# directly, the set_scaling_model function seems to run ok. 
# Maybe it is because the help files for set_scaling_model in 
# wrapper_functions.R have not been written yet.

set_scaling_model <- function(f0 = 0.6,
                              alpha = 0.4,
                              r_pp = 1e-1,
                              n = 2/3,
                              q =3/4,
                              kappa = 7e10,
                              erepro = 0.1,
                              beta = 100,
                              sigma = 1.3,
                              h = 30,
                              ks = 4,
                              no_sp = 11,
                              min_egg = 10^(-4),
                              min_w_pp = min_egg/(beta*exp(5*sigma)),
                              max_w = 10^3,
                              min_w_inf = 10,
                              min_w_mat = 10^(0.4),
                              no_w = log10(max_w/min_egg)*100+1){
  # Set exponents
  p <- n
  lambda <- 2+q-n
  # Set grid points and characteristic sizes 
  
  # The characteristic weights of the different species are defined by 
  # min_egg,min_w_mat, min_w_inf, max_w_inf and no_sp, in the sense 
  # that the egg weights of the no_sp species are logarithmically 
  # evenly spaced, ranging from min_w=min_egg to max_w=max_w_inf. 
  # The maturity weights of the species can be obtained by muliplying 
  # the egg_weights by min_w_mat/min_egg. The asymptotic weights 
  # of the species can be obtained by mulitiplying the egg weights by 
  # min_w_inf/min_egg.
  # The no_w weights, which we keep track of the abundance of fish at,
  # form a logarithmically evenly spaced vector w, ranging from 
  # min_w to max_w. The vector w_full is obtained by extending w 
  # down to the size range min_w_pp, where min_w_pp is chosen by default
  # to be so small that almost no fish can eat it.
  min_w <- min_egg
  # min_egg and max_w already lie on grid points in w. Let us round   
  # min_w_mat up to the nearest grid point.
  delt <- (log10(max_w)-log10(min_w))/(no_w-1)
  v <- min_w_mat
  j <- 1+ceiling((log10(v)-log10(min_w))/delt)
  v <- 10^(log10(min_w)+(j-1)*delt)
  min_w_mat <- v
  # Let us round min_w_inf up to the nearest grid point.
  v <- min_w_inf
  j <- 1+ceiling((log10(v)-log10(min_w))/delt)
  v <- 10^(log10(min_w)+(j-1)*delt)
  min_w_inf <- v
  # Determine maximum egg size
  max_egg <- max_w*min_egg/min_w_inf 
  log10_minimum_egg <- log10(min_egg)
  log10_maximum_egg <- log10(max_egg)
  # Determine logrithmic spacing of egg weights
  dist_sp <- (log10_maximum_egg-log10_minimum_egg)/(no_sp-1) 
  species <- 1:no_sp
  # Determine egg weights w_min for all species
  x_min <- seq(log10_minimum_egg, by = dist_sp, length.out = no_sp)
  w_min <- 10^x_min
  # Use ratios to determine w_inf and w_mat from w_min
  w_inf <- w_min*min_w_inf/min_egg
  w_mat <- w_min*min_w_mat/min_egg
  # Build Params Object 
  species_params <- data.frame(
    species = 1:no_sp,
    w_min = w_min,
    w_inf = w_inf,
    w_mat = w_mat,
    h = h,
    ks = ks,
    beta = beta,
    sigma = sigma,
    z0 = 0,
    alpha = alpha,
    erepro = erepro,
    sel_func = "knife_edge", # not used but required
    knife_edge_size = 1000
  )
  params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                        kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
                        min_w_pp = min_w_pp, w_pp_cutoff = max_w, r_pp = r_pp)
  # gamma is determined by mizerparams
  gamma <- params@species_params$gamma[1]
  w <- params@w
  # Get constants for analytic solution
  mu0 <- (1-f0) * sqrt(2*pi) * kappa * gamma * sigma *
    (beta^(n-1)) * exp(sigma^2 * (n-1)^2 / 2)
  hbar <- alpha * h * f0 - ks
  pow <- mu0/hbar/(1-n)
  n_mult <- (1 - (w/w_inf[1])^(1-n))^(pow-1) * (1 - (w_mat[1]/w_inf[1])^(1-n))^(-pow)
  n_mult[w < w_mat[1]] <- 1
  n_mult[w >= w_inf[1]] <- 0
  # Create steady state solution n_exact for species 1
  n_exact <- w  # Just to get array with correct dimensions and names
  n_exact <- ((w_min[1]/w)^(mu0/hbar) / (hbar * w^n)) * n_mult
  n_exact[w < w_min[1]] <- 0
  # Use n_exact as a template to create solution initial_n for all species
  initial_n <- params@psi
  initial_n[,] <- 0
  w_inf_idx <- w_inf
  for (i in 1:no_sp) {
    w_inf_idx[i] <- length(w[w<=w_inf[i]])
    initial_n[i, params@species_params$w_min_idx[i]:
                (params@species_params$w_min_idx[i]+
                   (w_inf_idx[1]-params@species_params$w_min_idx[1]))] <-
      n_exact[params@species_params$w_min_idx[1]:
                (params@species_params$w_min_idx[1]+
                   (w_inf_idx[1]-params@species_params$w_min_idx[1]))] *
      (w_min[1]/w_min[i])^lambda
  }
  # rescale fish abundance to line up with background resource spectrum
  v <- sqrt(min(w_mat)*max(w_mat))
  v_idx <- length(w[w<v])
  # The resulting steady state is n_output 
  n_output <- initial_n*(kappa*w[v_idx]^(-lambda))/sum(initial_n[,v_idx])
  # Setup plankton
  plankton_vec <- (kappa*w^(-lambda))-colSums(n_output)
  plankton_vec[plankton_vec<0] <- 0
  plankton_vec[min(which(plankton_vec==0)):length(plankton_vec)] <- 0
  params@cc_pp[sum(params@w_full<=w[1]):length(params@cc_pp)] <- plankton_vec
  initial_n_pp <- params@cc_pp
  # The cc_pp factor needs to be higher than the desired steady state in
  # order to compensate for predation mortality
  m2_background <- getM2Background(params, n_output, initial_n_pp)
  params@cc_pp <- (1+m2_background/params@rr_pp) * initial_n_pp
  # Setup background death and steplike psi
  m2 <- getM2(params, n_output, initial_n_pp)
  for (i in 1:no_sp) {
    params@psi[i, ] <- (w/w_inf[i])^(1-n)
    params@psi[i, w < (w_mat[i]-1e-10)] <- 0
    params@psi[i, w > (w_inf[i]-1e-10)] <- 1
    params@mu_b[i, ] <- mu0 * w^(n-1) - m2[i,]
  }
  # Set erepro to meet boundary condition
  rdi <- getRDI(params, n_output, initial_n_pp)
  gg <- getEGrowth(params, n_output, initial_n_pp)
  effort <- 0
  mumu <- getZ(params, n_output, initial_n_pp, effort = effort)
  erepro_final <- rdi
  for (i in (1:no_sp)){
    gg0 <- gg[i,params@species_params$w_min_idx[i]]
    mumu0 <- mumu[i,params@species_params$w_min_idx[i]]
    DW <- params@dw[params@species_params$w_min_idx[i]]
    erepro_final[i] <- erepro*(n_output[i,params@species_params$w_min_idx[i]]*(gg0+DW*mumu0))/rdi[i]
  }
  params@species_params$erepro <- erepro_final
  # Turn off R_max
  params@srr <- function(rdi, species_params) {return(rdi)}
  # Record abundance of fish and resources at steady state, as slots.
  params@initial_n <- n_output
  params@initial_n_pp <- initial_n_pp
  return(params)
}

params2 <- set_scaling_model()

t_max <- 5
sim <- project(params2, t_max=t_max, dt=0.01, t_save=t_max/100 ,effort = 0, 
               initial_n = params2@initial_n, initial_n_pp = params2@initial_n_pp)
plot(sim)