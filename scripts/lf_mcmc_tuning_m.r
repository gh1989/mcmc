library(zoo)
#setwd(sprintf("C:/cygwin64/home/Gregg/mcmc/output/%s",subfolder))
#par(mfrow=c(5,4) ,oma = c(0, 0, 2, 0) )
#par(mar=c(4,2,2,4))

###################################
# EXTRA DATA RATIO  PLOTS         #
###################################

tune_variable = "extra_data_ratio"
setwd(sprintf("C:/cygwin64/home/Gregg/mcmc/output/%s",tune_variable))

###################################
# Standard deviation              #
###################################

tune_variable = "extra_data_ratio"
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
  
  real_ess[i] = ess(real_parts)/100000
  imag_ess[i] = ess(imag_parts)/100000
  
  var_x[i] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
}

plot(var_x, real_var, log='xy', xlab='Data Ratio', ylab="Posterior Variance", main="", pch=2, col="red")
points(var_x, imag_var, log='xy', main="Im", pch=5, col="blue")

legend('topright',c("Re","Im"), lty=1, col=c('red', 'blue'), bty='n', cex=.75, pch=c(2,5))

###################################
# ESS                             #
###################################

plot(var_x, real_ess, log='xy', xlab='Data Ratio', ylab="ESS", main="", pch=2, col="red")
points(var_x, imag_ess, log='xy', main="Im", pch=5, col="blue")

legend('topright',c("Re","Im"), lty=1, col=c('red', 'blue'), bty='n', cex=.75, pch=c(2,5))

###################################
# Real part posterior for given K #
###################################
first_plot = 1
density_plots = c(2,3)
ptys=c(2,3)
num_k = length(density_plots) + 1

vals = c( first_plot, density_plots )

mcmc_data = read.table(list_of_files[first_plot])
x <- head( readLines(list_of_files[first_plot]))
meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
vals[1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
print(vals[1])
real_parts = mcmc_data[7:7]
plot(density(ts(real_parts)), pty=1, main="Posterior Re", xlim=c(0,0.55), ylim=c(0,25) )

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
legend('topright',as.character(vals), bty='n', 
       cex=.75, lty=c(1, 2, 3) )

###################################
# Imag part posterior for given K #
###################################
first_plot = 1
density_plots = c(2,3)
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
legend('topright',as.character(vals), bty='n', 
       cex=.75, lty=c(1, 2, 3) )


#mtext("Path Length", outer=TRUE, cex=1.5)

