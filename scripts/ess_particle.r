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
  x <- readLines(list_of_files[i], n=8)
  meta_data_string <- grep( pattern=the_pattern, x, value=TRUE)
  Ps[i] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  print(Ps[i])
  ess_fourth[i] = ess( mcmc_data[,7] ) + ess( mcmc_data[,8] )
  ess_fourth[ is.nan(ess_fourth) ] <- 1
  print(ess_fourth[i])
}

legend_vals =  c( expression(Q == 10^{1}),
                  expression(Q == 10^{2}),
                  expression(Q == 10^{3}))


colours = c("green","blue","red")

the_xlim <- c( min(Ps), max(Ps) )
the_ylim <- c( min(ess_fourth), max(ess_fourth) )

smoothingSpline = smooth.spline(ess_fourth[1:size_runs] ~ Ps[1:size_runs], spar=0.4, tol = 0.1)
plot( ess_fourth[1:size_runs] ~ Ps[1:size_runs], log="xy", xlab=the_var, ylab="Particles", col="green", pch=17, xlim=the_xlim, ylim=the_ylim)
lines(smoothingSpline, col="green")

smoothingSpline = smooth.spline(ess_fourth[(size_runs+1):(2*size_runs)] ~ Ps[(size_runs+1):(2*size_runs)], spar=0.1)
points( ess_fourth[(size_runs+1):(2*size_runs)] ~ Ps[(size_runs+1):(2*size_runs)], col="blue", pch=18 )
lines(smoothingSpline, col="blue")

smoothingSpline = smooth.spline(ess_fourth[(2*size_runs+1):(3*size_runs)] ~ Ps[(2*size_runs+1):(3*size_runs)] , spar=1.0)
points( ess_fourth[(2*size_runs+1):(3*size_runs)] ~ Ps[(2*size_runs+1):(3*size_runs)] , col="red", pch=19 )
lines(smoothingSpline, col="red")

legend( 'topright',
        legend_vals,
        bty="n",
        col=colours, 
        pch=c(17, 18, 19),
        lty=1,
        inset = c(0.1,0))
