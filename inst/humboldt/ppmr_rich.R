
library(tidyverse)
ppmr <- read_csv("inst/humboldt/Data_PPMR.csv")

a_ppmr <- filter(ppmr, Species == "anchovy")

log(a_ppmr$wprey)

data <- a_ppmr
low <- log(6.18e-9)
high <- log(1.68e-8) 


sum(data[(low<log(data[,2]))&(log(data[,2])<high),4])

#sum(data$Nprey[data[(low<log(data[,2]))&(log(data[,2])<high)])

add_up <- function(data,low,high){
    return(sum(data[(low<=log(data[,2]))&(log(data[,2])<high),4]))
}

no_bins <- 50
LEAST <- min(log(data$wprey))
MOST  <- max(log(data$wprey))

bins <- seq(from = min(log(data$wprey)), to = max(log(data$wprey)),length.out = no_bins)
bin_count <- bins
long_bins <- c(bins,bins[length(bins)]+(bins[length(bins)]-bins[length(bins)-1]))
for (t in (1:length(bins))){
    low <- long_bins[t]
    high <- long_bins[t+1]
    bin_count[t] <- sum(data$Nprey[(low<=log(data$wprey))&(log(data$wprey)<high)])
}

barplot(bin_count, width=bins)
bin_count
barplot(bin_count)
low
high
sum(data$Nprey[(low<=log(data$wprey))&(log(data$wprey)<high)])

barplot(bin_count *exp(bins))

