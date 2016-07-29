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
    
    plot(ts(pmcmcdata[2:2]))
    plot(rollmean(ts(pmcmcdata[2:2]),100))
    acf(pmcmcdata[2:2])
    hist(ts(pmcmcdata[2:2]), 30, main="")
    
}
library(zoo)
par(mar=c(4,2,2,4))
setwd("~/inferring_geometries/output")
data1 = read.table("mcmc_old_data_ts.txt")
pmcmcSummary2d(data1)
