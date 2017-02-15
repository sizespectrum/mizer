
  for (j in 1:noSpecies){
    Beta <- log(object@species_params$beta)[j]
    sigma <- object@species_params$sigma[j]
    
    
    
    Delta <- dx*round(min(2*sigma, Beta)/dx)
    Beta <- dx*round(Beta/dx)
    Delta <- Beta
    min_cannibal <- 1+floor((Beta-Delta)/dx)
    #min_cannibal <- 0
    
    P <- x[length(x)] + 2*Delta
    no_P <- 1+ceiling(P/dx)  # P/dx should already be integer 
    x_P <- (1:no_P)*dx#+Beta-Delta-dx
    
    phi <- rep(0, length(x_P))
    
    phi[abs(x_P+Beta-P)<Delta] <- exp(-(x_P[abs(x_P+Beta-P)<Delta] + Beta - P)^2/(2*sigma^2))
    
################################
    
    noSpecies <- dim(res@interaction)[1]
    no_Pvec <- rep(0, noSpecies)
    Pvec <- rep(0, noSpecies)
    for (j in 1:noSpecies){
      Beta <- log(res@species_params$beta)[j]
      sigma <- res@species_params$sigma[j]
      Delta <- dx*round(min(2*sigma, Beta)/dx)
      Beta <- dx*round(Beta/dx)
      Delta <- Beta
      min_cannibal <- 1+floor((Beta-Delta)/dx)
      P <- x[length(x)] + 2*Delta
      no_P <- 1+ceiling(P/dx)  # P/dx should already be integer 
      x_P <- (1:no_P)*dx#+Beta-Delta-dx
      no_Pvec[j] <- no_P
      Pvec[j] <- P
    }
    no_P <- max(no_Pvec)
    P <- max(Pvec)
    # # # 
    
    phiMortality <- matrix(0,nrow = noSpecies, ncol = no_P)
    for (j in 1:noSpecies){
      Beta <- log(res@species_params$beta)[j]
      sigma <- res@species_params$sigma[j]
      Delta <- dx*round(min(2*sigma, Beta)/dx)
      Beta <- dx*round(Beta/dx)
      Delta <- Beta
      min_cannibal <- 1+floor((Beta-Delta)/dx)
    #  P <- x[length(x)] + 2*Delta
    #  no_P <- 1+ceiling(P/dx)  # P/dx should already be integer 
      x_P <- (1:no_P)*dx#+Beta-Delta-dx
      phi <- rep(0, length(x_P))
      phi[abs(x_P+Beta-P)<Delta] <- exp(-(x_P[abs(x_P+Beta-P)<Delta] + Beta - P)^2/(2*sigma^2)) 
      phiMortality[j, 1:length(phi)] <- phi
    }