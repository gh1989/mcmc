# Outputs for pMCMC time series.
pmcmcSummary <- function( pmcmcdata )
{
    par(mfrow = c(3, 1))
    plot(ts(pmcmcdata))
    acf(pmcmcdata)
    hist(ts(pmcmcdata), 30)
    
}

pmcmcSummary2d <- function( pmcmcdata )
{
    par(mfrow = c(4, 1))
    
    plot(ts(pmcmcdata))
    plot(rollmean(ts(pmcmcdata),50))
    acf(pmcmcdata)
    hist(ts(pmcmcdata), 30, main="")
    
}
library(zoo)
par(mar=c(4,2,2,4))
setwd("~/mcmc/output")
pmcmcdata_timeseries = read.table("ou_mcmc_timeseries.txt")
pmcmcSummary2d(pmcmcdata_timeseries)
