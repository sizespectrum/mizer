library(tidyverse)
ppmr <- read_csv("Data_PPMR.csv")
ppmr
lambda <- 2

moments <- ppmr %>%
    group_by(Species) %>%
    summarise(totalNprey = sum(Nprey),
              m1 = sum(log(wprey) * Nprey / totalNprey), 
              m2 = sum(log(wprey)^2 * Nprey / totalNprey),
              var = m2 - m1^2,
              mu = m1 + var * (lambda - 1),
              sigma = sqrt(var),
              logbeta = log(first(wpredator)) - mu,
              beta = exp(logbeta)
    )
View(moments)


dens <- function(x, mu, sigma, lambda) {
    exp(-x * (lambda - 1) - (x - mu)^2 / (2 * sigma^2)) / 
        (sqrt(2 * pi) * sigma) / 
        exp((lambda - 1) * (-mu + (lambda - 1) * sigma^2 / 2))
}

s_ppmr <- filter(ppmr, Species == "sardine")
exploded <- data.frame(logwprey = rep(log(s_ppmr$wprey), 
                                      round(s_ppmr$Nprey * 1000)))

ggplot(exploded) +
    geom_histogram(aes(logwprey, stat(density)), bins = 20) +
    stat_function(fun = dens, args = list(mu = moments$mu[2],
                                          sigma = moments$sigma[2],
                                          lambda = lambda)
    )
