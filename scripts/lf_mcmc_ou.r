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
    ou_ts <- ts(mcmc_data[1:1])
    plot(ou_ts)
    plot(rollmean(ou_ts,50))
    acf(ou_ts)
    hist(ou_ts, 30, main="")
    
}
plotHistogramFromFilenameOU<-function(f)
{
  ou_mcmc_data = read.table(f)
  mcmcSummary2d(ou_mcmc_data)
}

library(zoo)
par(mar=c(4,2,2,4))
#setwd("C:/cygwin64/home/Gregg/mcmc/output")
list_of_files = list.files(pattern=".*TimeSeries_.*.txt")
lapply(list_of_files, plotHistogramFromFilenameOU)

