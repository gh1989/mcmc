library(mcmcse)

setwd("C:/cygwin64/home/Gregg/mcmc/output/PMCMC")
setwd(folder)
setwd(the_var)

list_of_files = list.files(pattern="pMCMCTimeSeries_.*.txt")
len = length(list_of_files)
Ps = numeric(len)
ess_fourth = numeric(len)

for(i in 1:len)
{
  mcmc_data = read.table(list_of_files[i])
  x <- readLines(list_of_files[i], n=9)
  meta_data_string <- grep( pattern="_number_particles", x, value=TRUE)
  N = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  
  
  #meta_data_string <- grep( pattern="_single_mode", x, value=TRUE)
  #single <- as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  #if (!single)
  #{
  #  Ps[i] = 2*M*(M+1)
  #}
  
  Ps[i] = N

  print(Ps[i])
  ess_fourth[i] = ess( mcmc_data[,7] ) + ess( mcmc_data[,8] )
  ess_fourth[ is.nan(ess_fourth) ] <- 1
  print(ess_fourth[i])
}

legend_vals =  c( expression(dim(C) == 1),
                  expression(dim(C) == 4),
                  expression(dim(C) == 12),
                  expression(dim(C) == 40) )


colours = c("green","blue","red", "purple")

the_xlim <- c( min(Ps), max(Ps) )
the_ylim <- c( min(ess_fourth), max(ess_fourth) )

smoothingSpline = smooth.spline(ess_fourth[1:size_runs] ~ Ps[1:size_runs], spar=0.4, tol = 0.1)
plot( ess_fourth[1:size_runs] ~ Ps[1:size_runs], log="xy", xlab="Q", ylab="Particles", col="green", pch=17, xlim=the_xlim, ylim=the_ylim)
lines(smoothingSpline, col="green")

smoothingSpline = smooth.spline(ess_fourth[(size_runs+1):(2*size_runs)] ~ Ps[(size_runs+1):(2*size_runs)], spar=0.1)
points( ess_fourth[(size_runs+1):(2*size_runs)] ~ Ps[(size_runs+1):(2*size_runs)], col="blue", pch=18 )
lines(smoothingSpline, col="blue")

smoothingSpline = smooth.spline(ess_fourth[(2*size_runs+1):(3*size_runs)] ~ Ps[(2*size_runs+1):(3*size_runs)] , spar=1.0)
points( ess_fourth[(2*size_runs+1):(3*size_runs)] ~ Ps[(2*size_runs+1):(3*size_runs)] , col="red", pch=19 )
lines(smoothingSpline, col="red")

smoothingSpline = smooth.spline(ess_fourth[(3*size_runs+1):(4*size_runs)] ~ Ps[(3*size_runs+1):(4*size_runs)] , spar=1.0)
points( ess_fourth[(3*size_runs+1):(4*size_runs)] ~ Ps[(3*size_runs+1):(4*size_runs)] , col="purple", pch=19 )
lines(smoothingSpline, col="purple")

legend( 'topright',
        legend_vals,
        bty="n",
        col=colours, 
        pch=c(17, 18, 19, 20),
        lty=1,
        inset = c(0.1,0))
