# Outputs for MCMC time series plots.
mcmcSummary <- function( mcmc_data )
{
    par(mfrow = c(3, 1))
    plot(ts(mcmc_data))
    acf(mcmc_data)
    hist(ts(mcmc_data), 30)
    
}

mcmcSummary2d <- function( mcmc_data )
{
    par(mfrow = c(4, 1))
    
    plot(ts(mcmc_data))
    plot(rollmean(ts(mcmc_data),50))
    acf(mcmc_data)
    hist(ts(mcmc_data), 30, main="")
    
}
library(zoo)
par(mar=c(4,2,2,4))
setwd("~/mcmc/output")
ou_mcmc_data = read.table("ou_mcmc_timeseries.txt")
mcmcSummary2d(ou_mcmc_data)
