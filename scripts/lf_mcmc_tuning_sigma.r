#################
# PLOTS         #
#################
algo_type = "LFMCMC"
setwd(sprintf("C:/cygwin64/home/Gregg/mcmc/output/%s/sigma_only/%s/",algo_type, tune_variable))

###################################
# Standard deviation              #
###################################

list_of_files = list.files(pattern=".*TimeSeries_.*.txt")
len = length(list_of_files)
real_var = numeric(len)
var_x = numeric(len)
real_ess = numeric(len)

for(i in 1:len)
{
  mcmc_data = read.table(list_of_files[i])
  x <- head( readLines(list_of_files[i]))
  meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
  real_parts = mcmc_data[9:9]
  real_var[i] = var(as.numeric(unlist(real_parts)))
  real_ess[i] = ess(real_parts)/mcmc_trials
  var_x[i] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  print( var_x[i] )
}

ylim <- c( min(real_var), max(real_var) )

plot(var_x, real_var, log=log_xy, xlab=nice_name, ylab="Posterior Variance", main="Posterior Var.", pch=2, col="red", ylim=ylim)
smoothingSpline = smooth.spline(real_var~var_x, spar=0.4)
lines(smoothingSpline, col="red")

###################################
# ESS                             #
###################################
ylim <- c( min(real_ess), max(real_ess) )
plot(var_x, real_ess, log=log_xy, xlab=nice_name, ylab="ESS", pch=2, col="red", ylim=ylim, main="ESS")

smoothingSpline = smooth.spline(real_ess~var_x, spar=0.4)
lines(smoothingSpline, col="red")

###################################
# Real part posterior for given K #
###################################
num_k = length(density_plots) + 1

vals = c( first_plot, density_plots )

mcmc_data = read.table(list_of_files[first_plot])
x <- head( readLines(list_of_files[first_plot]))
meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
vals[1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
print(vals[1])
real_parts = mcmc_data[9:9]
plot(density(ts(real_parts)), pty=1, main="Posterior", xlim=c(-8.0, -5.0) )
abline(v=-6.9, lty=12, col="red")

for(i in 1:length(density_plots) )
{  
  idx = density_plots[i]
  mcmc_data = read.table(list_of_files[idx])
  real_parts = mcmc_data[9:9]
  x <- head( readLines(list_of_files[idx]))
  meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
  lines(density(ts(real_parts), adjust = 9), lty=i+1)
  vals[i+1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  print(vals[i+1])
}
legend('topright',as.character(vals), bty='n', xpd=TRUE, cex=.75, lty=c(1, 2, 3) )

