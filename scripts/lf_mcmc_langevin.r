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
    hist(ts(mcmc_data[7:7]), 60, main="Re(mode 4)")
    hist(ts(mcmc_data[8:8]), 60, main="Re(mode 4)")
}

histgramsModes<-function( mcmc_data )
{
    par(mfrow = c(4, 2))
    for(i in 1:4)
    {
        hist(ts(mcmc_data[(i*2-1):(i*2-1)]), 60, main=sprintf("Re(mode %d)",i))      
        hist(ts(mcmc_data[(i*2):(i*2)]), 60, main=sprintf("Im(mode %d)",i))      
    }
}

sigmaSummary <- function( mcmc_data )
{
  par(mfrow = c(4, 1))
  plot(ts(mcmc_data[9:9]))
  plot(rollmean(ts(mcmc_data[9:9]),100))
  acf(mcmc_data[9:9])
  print(mean(ts(mcmc_data[9:9])))
  hist( ts(mcmc_data[9:9]), breaks=25, freq=FALSE, main="Posterior Distribution")
}


plotHistogramFromFilename<-function(f)
{
  mcmc_data = read.table(f)
  histgramsModes(mcmc_data)
  mcmcSummary2d(mcmc_data)
  #sigmaSummary(mcmc_data)
}
  
library(zoo)
par(mar=c(4,2,2,4))
setwd("C:/cygwin64/home/Gregg/mcmc/output")
list_of_files = list.files(pattern="LangevinTimeSeries_.*.txt")
lapply(list_of_files, plotHistogramFromFilename)
