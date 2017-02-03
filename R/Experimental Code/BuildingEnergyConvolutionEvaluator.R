



#####################

install.packages("mizer")

library(mizer)

params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7)

params@interaction

summary(params)

params@pred_kernel

#####################

#object@pred_kernel[i,w,wp]=Phi_i(wp/w)

# perhaps we can use the above information directly, anyway we can use

wFull <- object@w_full 

w0 <- eggsize

# how to find eggsize ?

xFull <- log(wFull/w0)

# how to find sigma and Beta ? note, our Beta = log(beta_from_vignette) ?

s <- exp(-(xFull -Beta)^2/(2*sigma^2))

N <- length(s)

dx <- xFull[2]-xFull[1]

# n[j,wp]= N_j(w_p)

f <- sweep(sweep((object@interaction %*% n), 1, n_pp, "+"), 2, wFull, "*")[i]

# could it be f <- sweep(sweep((object@interaction %*% n), 1, n_pp, "+"), 2, wFull^2, "*")[i] ?

# Alternative loop based approach...

A <- (object@interaction %*% n)
for (i in 1:(dim(object@interaction)[1])){
 A[i,] <- A[i,] + nn_p}

for (i in 1:(dim(object@interaction)[2])){
    A[,i] <- A[,i] * wFull}

f <- A[i, ]


# can/should we do the convolution integral for all species i simultaineously ?

energy <- dx*Re(fft(fft(s)*fft(f), inverse=TRUE)/N)

