#' ---
#' title: "A Stable Community"
#' author: "Richard Southwell"
#' output: pdf_document
#' ---


library(mizer)

#' We use the analytic solution to our trait based model, with variable egg size, to set up a 
#' community thats aggregation produces a background power law size spectrum that is self stableizing.
#' The stablity of this fixed point is promoted by our inclusion of density dependence.
#' 

#' ## Without density dependence
# ----
#' ### Set parameters 
#' Global Parameters
#' 

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
                              #dist_sp = 0.2,
                              no_sp = 11,
                              log10_minimum_egg = -4,
                              log10_maximum_egg = -2,
                              min_w_pp = 1e-8){
  
  
  p <- n
  lambda <- 2+q-n
  
  #' Species Parameters
  
  # ----
  #' ### Set grid points and characteristic sizes 
  
  #! we should describe size in terms of min_egg,min_w_mat,min_w_inf,
  # max_w_inf,no_sp and no_sp, so the 5, 4.5, dist_sp,log10_minimum_egg,
  # log10_maximum_egg need to be replaced.
  
  #! how much freedom do we have over the choice of no_w ?
  
  #! set min_w_pp properly later, by considering beta, sigma,min_w
  
  ##!! altered things so that no_sp is defined, and dist_sp is 
  # computed as a result. Notice new code has issues if no_sp=1
  ## no_sp <- (log10_maximum_egg-log10_minimum_egg)/dist_sp + 1
  dist_sp <- (log10_maximum_egg-log10_minimum_egg)/(no_sp-1) 
  
  
  species <- 1:no_sp
  x_min <- seq(log10_minimum_egg, by = dist_sp, length.out = no_sp)
  w_min <- 10^x_min
  w_inf <- 10^(x_min+5)
  w_mat <- 10^(x_min+4.4)  # This is about a quarter of w_inf
  min_w <- min(w_min)
  max_w <- max(w_inf)
  no_w <- log10(max_w/min_w)*100+1
  
  # ----
  #' ### Build Params Object 
  
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
  
  #params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
  #                      kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
  #                      min_w_pp = min_w_pp, w_pp_cutoff = max_w, r_pp = r_pp,
  #                      chi = 0)
  params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                        kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
                        min_w_pp = min_w_pp, w_pp_cutoff = max_w, r_pp = r_pp)
  
  #' gamma is determined by mizerparams. Note that density dependence is currently off
  gamma <- params@species_params$gamma[1]
  w <- params@w
  
  
  # ----
  #' ### Determine analytic solution
  
  
  mu0 <- (1-f0) * sqrt(2*pi) * kappa * gamma * sigma *
    (beta^(n-1)) * exp(sigma^2 * (n-1)^2 / 2)
  hbar <- alpha * h * f0 - ks
  pow <- mu0/hbar/(1-n)
  n_mult <- (1 - (w/w_inf[1])^(1-n))^(pow-1) * (1 - (w_mat[1]/w_inf[1])^(1-n))^(-pow)
  n_mult[w < w_mat[1]] <- 1
  n_mult[w >= w_inf[1]] <- 0
  n_exact <- w  # Just to get array with correct dimensions and names
  n_exact <- ((w_min[1]/w)^(mu0/hbar) / (hbar * w^n)) * n_mult
  n_exact[w < w_min[1]] <- 0
  
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
  
  v <- sqrt(min(w_mat)*max(w_mat))
  v_idx <- length(w[w<v])
  #n_output <- initial_n*(kappa*w[v_idx]^(-lambda))/sum(sqrt(initial_n[,v_idx]*initial_n[,v_idx+19]))
  n_output <- initial_n*(kappa*w[v_idx]^(-lambda))/sum(initial_n[,v_idx])
  
  # more generally, just sum over 1:no_background species
  
  
  # ----
  #' ### Setup plankton
  
  plankton_vec <- (kappa*w^(-lambda))-colSums(n_output)
  plankton_vec[plankton_vec<0] <- 0
  plankton_vec[min(which(plankton_vec==0)):length(plankton_vec)] <- 0
  params@cc_pp[sum(params@w_full<=w[1]):length(params@cc_pp)] <- plankton_vec
  initial_n_pp <- params@cc_pp
  # The cc_pp factor needs to be higher than the desired steady state in
  # order to compensate for predation mortality
  m2_background <- getM2Background(params, n_output, initial_n_pp)
  params@cc_pp <- (1+m2_background/params@rr_pp) * initial_n_pp
  
  # ----
  #' ### Setup background death and steplike psi
  
  m2 <- getM2(params, n_output, initial_n_pp)
  
  for (i in 1:no_sp) {
    params@psi[i, ] <- (w/w_inf[i])^(1-n)
    params@psi[i, w < (w_mat[i]-1e-10)] <- 0
    params@psi[i, w > (w_inf[i]-1e-10)] <- 1
    params@mu_b[i, ] <- mu0 * w^(n-1) - m2[i,]
  }
  
  # ----
  #' ### Set erepro to meet boundary condition
  
  rdi <- getRDI(params, n_output, initial_n_pp)
  gg <- getEGrowth(params, n_output, initial_n_pp)
  effort <- 0
  mumu <- getZ(params, n_output, initial_n_pp, effort = effort)
  erepro_final <- rdi
  for (i in (1:no_sp)){
    #  erepro_final[i] <- erepro*(gg[i,params@species_params$w_min_idx[i]]*n_output[i,params@species_params$w_min_idx[i]])/
    #    rdi[i]
    gg0 <- gg[i,params@species_params$w_min_idx[i]]
    mumu0 <- mumu[i,params@species_params$w_min_idx[i]]
    DW <- params@dw[params@species_params$w_min_idx[i]]
    erepro_final[i] <- erepro*(n_output[i,params@species_params$w_min_idx[i]]*(gg0+DW*mumu0))/rdi[i]
  }
  
  params@species_params$erepro <- erepro_final
  
  # ----
  #' ### Simulate without Rmax or chi
  
  params@srr <- function(rdi, species_params) {return(rdi)}
  #params@chi <- 0.0
  
  params@initial_n <- n_output
  params@initial_n_pp <- initial_n_pp
  
  
  return(params)
}

params2 <- set_scaling_model()

t_max <- 5
sim <- project(params2, t_max=t_max, dt=0.01, t_save=t_max/100 ,effort = 0, 
               initial_n = params2@initial_n, initial_n_pp = params2@initial_n_pp)
plot(sim)