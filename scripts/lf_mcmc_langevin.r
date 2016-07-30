# Outputs for MCMC time series plots.
mcmcSummary2d <- function( mcmc_data )
{
    par(mfrow = c(4, 2))
    plot(ts(mcmc_data[7:7]))
    plot(ts(mcmc_data[8:8]))
    plot(rollmean(ts(mcmc_data[7:7]),100))
    plot(rollmean(ts(mcmc_data[8:8]),100))
    acf(mcmc_data[7:7])
    acf(mcmc_data[8:8])
    hist(ts(mcmc_data[7:7]), 60, main="")
    hist(ts(mcmc_data[8:8]), 60, main="")
}

histgramsModes<-function( mcmc_data )
{
    par(mfrow = c(4, 2))
    for(i in 1:4)
    {
        hist(ts(mcmc_data[(i*2-1):(i*2-1)]), 60, main="")
        hist(ts(mcmc_data[(i*2):(i*2)]), 60, main="")      
    }
}

library(zoo)
par(mar=c(4,2,2,4))
setwd("~/mcmc/output")
mcmc_data = read.table("langevin_mcmc_timeseries.txt")
histgramsModes(mcmc_data)
