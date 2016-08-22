#################
# PLOTS         #
#################
mcmc_trials = 100000
setwd(sprintf("C:/cygwin64/home/Gregg/mcmc/output/LFMCMC/%s",tune_variable))

###################################
# Standard deviation              #
###################################

list_of_files = list.files(pattern=".*TimeSeries_.*.txt")
len = length(list_of_files)
real_var = numeric(len)
imag_var = numeric(len)
var_x = numeric(len)
real_ess = numeric(len)
imag_ess = numeric(len)

for(i in 1:len)
{
  mcmc_data = read.table(list_of_files[i])
  x <- head( readLines(list_of_files[i]))
  meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
  real_parts = mcmc_data[7:7]
  imag_parts = mcmc_data[8:8]
  real_var[i] = var(as.numeric(unlist(real_parts)))
  imag_var[i] = var(as.numeric(unlist(imag_parts)))
  
  real_ess[i] = ess(real_parts)/mcmc_trials
  imag_ess[i] = ess(imag_parts)/mcmc_trials
  
  var_x[i] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  print( var_x[i] )
}
ylim <- c( min(real_var,imag_var), max(real_var, imag_var) )
plot(var_x, real_var, log=log_xy, xlab=nice_name, ylab="Posterior Variance", main="", pch=2, col="red", ylim=ylim)
points(var_x, imag_var, log=log_xy, main="Im", pch=5, col="blue")
legend('topright',c("Re","Im"), lty=1, xpd=TRUE, inset=c(-0.25,0), col=c('red', 'blue'), bty='n', cex=.75, pch=c(2,5))

###################################
# ESS                             #
###################################
ylim <- c( min(real_ess,imag_ess), max(real_ess, imag_ess) )
plot(var_x, real_ess, log=log_xy, xlab=nice_name, ylab="ESS", main="", pch=2, col="red", ylim=ylim)
points(var_x, imag_ess, log=log_xy, main="Im", pch=5, col="blue")

legend('topright',c("Re","Im"), lty=1, xpd=TRUE, inset=c(-0.25,0), col=c('red', 'blue'), bty='n', cex=.75, pch=c(2,5))

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
real_parts = mcmc_data[7:7]
plot(density(ts(real_parts)), pty=1, main="Posterior Re", xlim=c(0,0.55), ylim=c(0,25), xpd=TRUE, inset=c(-0.25,0) )

for(i in 1:length(density_plots) )
{  
  idx = density_plots[i]
  mcmc_data = read.table(list_of_files[idx])
  real_parts = mcmc_data[7:7]
  x <- head( readLines(list_of_files[idx]))
  meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
  lines(density(ts(real_parts), adjust = 9), lty=i+1)
  vals[i+1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  print(vals[i+1])
}
legend('topright',as.character(vals), bty='n', xpd=TRUE, inset=c(-0.30,0),
       cex=.75, lty=c(1, 2, 3) )

###################################
# Imag part posterior for given K #
###################################
num_k = length(density_plots) + 1

vals = c( first_plot, density_plots )

mcmc_data = read.table(list_of_files[first_plot])
x <- head( readLines(list_of_files[first_plot]))
meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
vals[1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
print(vals[1])
imag_parts = mcmc_data[8:8]
plot(density(ts(imag_parts)), pty=1, main="Posterior Im", xlim=c(-0.55,0), ylim=c(0,25) )

for(i in 1:length(density_plots) )
{  
  idx = density_plots[i]
  mcmc_data = read.table(list_of_files[idx])
  imag_parts = mcmc_data[8:8]
  x <- head( readLines(list_of_files[idx]))
  meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
  lines(density(ts(imag_parts), adjust = 8), lty=i+1)
  vals[i+1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  print(vals[i+1])
}
legend('topright',as.character(vals), bty='n', cex=.75, lty=c(1, 2, 3), xpd=TRUE, inset=c(-0.3,0))


#mtext("Path Length", outer=TRUE, cex=1.5)

