# Outputs for MCMC time series plots.
library(zoo)

#setwd("C:/cygwin64/home/Gregg/mcmc/output/parallel_paths/1 param")
#setwd("~/mcmc/output")

mcmcSummary2d <- function( mcmc_data )
{
    par(mfrow = c(4, 2))
    plot(ts(mcmc_data[7:7]), xlab="Iteration", ylab="")
    plot(ts(mcmc_data[8:8]), xlab="Iteration", ylab="")
    plot(rollmean(ts(mcmc_data[7:7]),20), xlab="Iteration", ylab="")
    plot(rollmean(ts(mcmc_data[8:8]),20), xlab="Iteration", ylab="")
    acf(mcmc_data[7:7], ylab="")
    acf(mcmc_data[8:8], ylab="")
    plot( density(ts(mcmc_data[7:7])), main="Posterior Distribution", xlab="Re(c)", xlim=c(0.35,0.60))
    plot( density(ts(mcmc_data[8:8])), main="Posterior Distribution", xlab="Im(c)", xlim=c(-0.35,-0.60))
    
    print(mean(ts(mcmc_data[7:7])))
    print(sd(ts(mcmc_data[7:7])))  
    
    print(mean(ts(mcmc_data[8:8])))
    print(sd(ts(mcmc_data[8:8]))) 
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
  #par(mfrow = c(4, 3))
  mcmc_ts <- ts(mcmc_data[9:9])
  #plot(mcmc_ts, xlab="Iteration", ylab="", main="Time Series")
  plot(rollmean(ts(mcmc_data[9:9]),100), main="Average")
  acf(mcmc_data[9:9], main="ACF")
  #print(mean(ts(mcmc_data[9:9])))
  plot( density(mcmc_ts, adjust = 4), main="Posterior Distribution", xlim=c(-10, -2))
  print(mean(ts(mcmc_data[9:9])))
  print(sd(ts(mcmc_data[9:9])))
}


plotHistogramFromFilename<-function(f)
{
  mcmc_data = read.table(f)
  #histgramsModes(mcmc_data)
  #mcmcSummary2d(mcmc_data)
  sigmaSummary(mcmc_data)
}
  
par(mfrow = c(3, 3))
par(mar=c(4,2,2,4))

list_of_files = list.files(pattern="LFTimeSeries_.*.txt")
lapply(list_of_files, plotHistogramFromFilename)
