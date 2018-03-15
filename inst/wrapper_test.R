params2 <- set_scaling_model()

t_max <- 5
sim <- project(params2, t_max=t_max, dt=0.01, t_save=t_max/100 ,effort = 0, 
               initial_n = params2@initial_n, initial_n_pp = params2@initial_n_pp)
plot(sim)