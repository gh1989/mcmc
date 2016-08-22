# Outputs for MCMC time series plots.
library(zoo)

#setwd("C:/cygwin64/home/Gregg/mcmc/output/parallel_paths/1 param")
#setwd("~/mcmc/output")

mcmcSummary2d <- function( mcmc_data )
{
    par(mfrow = c(4, 2))
    plot(ts(mcmc_data[7:7]))
    plot(ts(mcmc_data[8:8]))
    plot(rollmean(ts(mcmc_data[7:7]),20))
    plot(rollmean(ts(mcmc_data[8:8]),20))
    acf(mcmc_data[7:7])
    acf(mcmc_data[8:8])
    hist(ts(mcmc_data[7:7]), 20, main="Re(mode 4)", xlab="Re(c)")#, xlim=c(0.30,0.70))
    hist(ts(mcmc_data[8:8]), 20, main="Im(mode 4)", xlab="Im(c)")#, xlim=c(-0.30,-0.70))
    
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
  par(mfrow = c(4, 1))
  plot(ts(mcmc_data[9:9]))
  plot(rollmean(ts(mcmc_data[9:9]),100))
  acf(mcmc_data[9:9])
  print(mean(ts(mcmc_data[9:9])))
  hist( ts(mcmc_data[9:9]), breaks=25, freq=FALSE, main="Posterior Distribution")
  #print(mean(ts(mcmc_data[9:9])))
  #print(sd(ts(mcmc_data[9:9])))
}


plotHistogramFromFilename<-function(f)
{
  mcmc_data = read.table(f)
  #histgramsModes(mcmc_data)
  mcmcSummary2d(mcmc_data)
  #sigmaSummary(mcmc_data)
}
  

par(mar=c(4,2,2,4))

list_of_files = list.files(pattern=".*TimeSeries_.*.txt")
lapply(list_of_files, plotHistogramFromFilename)
