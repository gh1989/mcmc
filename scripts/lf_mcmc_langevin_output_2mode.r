# Outputs for MCMC time series plots.
library(zoo)

#setwd("C:/cygwin64/home/Gregg/mcmc/output/parallel_paths/1 param")
#setwd("~/mcmc/output")

mcmcSummary2d <- function( f )
{
  
    mcmc_data = read.table(f)
    x <- readLines(f,9)
    meta_data_string <- grep( pattern="real_sigma", x, value=TRUE)
    real_sigma = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  
    print(real_sigma)
    
    par(mfrow = c(3, 2))
    #plot(ts(mcmc_data[7:7]), xlab="Iteration", ylab="")
    #plot(ts(mcmc_data[8:8]), xlab="Iteration", ylab="")
    plot(rollmean(ts(mcmc_data[7:7]),20), xlab="Iteration", ylab="", ylim=c(-1,1), main="Average Re(c)")
    abline(h=0.5, lty=12, col="red")
    plot(rollmean(ts(mcmc_data[8:8]),20), xlab="Iteration", ylab="", ylim=c(-1,1), main="Average Im(c)")
    abline(h=-0.5, lty=12, col="red")
    #plot(rollmean(ts(mcmc_data[9:9]),20), xlab="Iteration", ylab="", main=expression(paste("Average ",log(sigma[d]))))
    #abline(h=log(real_sigma), lty=12, col="red")
    
    acf(mcmc_data[7:7], main="")
    acf(mcmc_data[8:8], main="")
    #acf(mcmc_data[9:9], main="")
    
    plot( density(ts(mcmc_data[7:7])), main="Posterior Distribution", xlab="Re(c)", xlim=c(-1,1))
    abline(v=0.5, lty=12, col="red")
    plot( density(ts(mcmc_data[8:8])), main="Posterior Distribution", xlab="Im(c)", xlim=c(-1,1))
    abline(v=-0.5, lty=12, col="red")
    #plot( density(ts(mcmc_data[9:9]), adjust = 4), main="Posterior Distribution", xlab=expression(log(sigma[d])))
    #abline(v=log(real_sigma), lty=12, col="red")
    
    print(mean(ts(mcmc_data[7:7])))
    print(sd(ts(mcmc_data[7:7])))  
    print(mean(ts(mcmc_data[8:8])))
    print(sd(ts(mcmc_data[8:8]))) 
    print(mean(ts(mcmc_data[9:9])))
    print(sd(ts(mcmc_data[9:9]))) 
}

histgramsModes<-function( f )
{
    par(mfrow = c(4, 4))
  
    mcmc_data = read.table(f)
    x <- readLines(f,9)
    meta_data_string <- grep( pattern="real_sigma", x, value=TRUE)
    real_sigma = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
    
    print(real_sigma)
  
    true_vals_real = c(0.5,0.25,1.0,-1.0)
    true_vals_imag = c(-0.5,0.5,1.0,-1.7)
    
    modes = c( 0, 1, 2, 3 )
    
    real_title = c( expression(Re(c[1][","][0])),expression(Re(c[-1][","][1])),expression(Re(c[0][","][1])),expression(Re(c[1][","][1])) )
    imag_title = c( expression(Im(c[1][","][0])),expression(Im(c[-1][","][1])),expression(Im(c[0][","][1])),expression(Im(c[1][","][1])) )
    
    for(i in 1:4)
    {
        mcmc_ts <- ts(mcmc_data[(i*2-1):(i*2-1)])
        plot( density(mcmc_ts, adjust = 4), main=real_title[i] )
        abline(v=true_vals_real[i], lty=12, col="red")
        acf( mcmc_ts, main="" )
        mcmc_ts <-ts(mcmc_data[(i*2):(i*2)])
        plot( density(mcmc_ts, adjust = 4), main=imag_title[i] )
        abline(v=true_vals_imag[i], lty=12, col="red")
        acf( mcmc_ts, main="" )
    }
}

sigmaSummary <- function( f )
{
  #par(mfrow = c(4, 3))
  mcmc_data = read.table(f)
  mcmc_ts <- ts(mcmc_data[9:9])
  #plot(mcmc_ts, xlab="Iteration", ylab="", main="Time Series")
  
  x <- readLines(f,9)
  meta_data_string <- grep( pattern="real_sigma", x, value=TRUE)
  real_sigma = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  
  plot(rollmean(ts(mcmc_data[9:9]),100), main="Average", xlim=c(1000, 10000),ylim=c(log(real_sigma)-2.0,log(real_sigma)+2.0) )
  abline(h=log(real_sigma), lty=12, col="red")
  acf(mcmc_data[9:9], main="")
  #print(mean(ts(mcmc_data[9:9])))
  plot( density(mcmc_ts, adjust = 4), main="Posterior Distribution", xlim=c(log(real_sigma)-2.0,log(real_sigma)+2.0) )
  abline(v=log(real_sigma), lty=12, col="red")
  print(mean(ts(mcmc_data[9:9])))
  print(sd(ts(mcmc_data[9:9])))
}


plotHistogramFromFilename<-function(f)
{
  mcmc_data = read.table(f)
  #histgramsModes(f)
  mcmcSummary2d(f)
  #sigmaSummary(f)
}
  
par(mfrow = c(3, 3))
par(mar=c(4,2,2,4))

list_of_files = list.files(pattern=".*TimeSeries_.*.txt")
lapply(list_of_files, plotHistogramFromFilename)
