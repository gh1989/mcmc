library(zoo)
#setwd(sprintf("C:/cygwin64/home/Gregg/mcmc/output/%s",subfolder))

###################################
# PATH LENGTH     PLOTS           #
###################################

tune_variable = "path_length"
setwd(sprintf("C:/cygwin64/home/Gregg/mcmc/output/%s",tune_variable))

###################################
# Standard deviation              #
###################################

tune_variable = "path_length"
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

plot(var_x, real_var, log='xy', xlab='Path Length', ylab="Posterior Variance", main="", pch=2, col="red")
points(var_x, imag_var, log='xy', main="Im", pch=5, col="blue")

legend('topright',c("Re","Im"), lty=1, col=c('red', 'blue'), bty='n', cex=.75, pch=c(2,5))

###################################
# ESS                             #
###################################

plot(var_x, real_ess, log='xy', xlab='Path Length', ylab="ESS", main="", pch=2, col="red")
points(var_x, imag_ess, log='xy', main="Im", pch=5, col="blue")

legend('topright',c("Re","Im"), lty=1, col=c('red', 'blue'), bty='n', cex=.75, pch=c(2,5))

###################################
# Real part posterior for given K #
###################################
first_plot = 9
density_plots = c(3,5)
num_k = length(density_plots) + 1

vals = numeric(num_k)

mcmc_data = read.table(list_of_files[first_plot])
x <- head( readLines(list_of_files[first_plot]))
meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
vals[1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
print(vals[1])
real_parts = mcmc_data[7:7]
plot(density(ts(real_parts)), pty=1, main="Posterior Re" )

for(i in 1:length(density_plots) )
{  
  idx = density_plots[i]
  mcmc_data = read.table(list_of_files[idx])
  real_parts = mcmc_data[7:7]
  x <- head( readLines(list_of_files[idx]))
  meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
  lines(density(ts(real_parts)), lty=i+1)
  vals[i+1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  print(vals[i+1])
}
legend('topright',as.character(vals), bty='n', 
       cex=.75, lty=c(1, density_plots) )

###################################
# Imag part posterior for given K #
###################################
first_plot = 9
density_plots = c(3,5)
num_k = length(density_plots) + 1

vals = numeric(num_k)

mcmc_data = read.table(list_of_files[first_plot])
x <- head( readLines(list_of_files[first_plot]))
meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
vals[1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
print(vals[1])
imag_parts = mcmc_data[8:8]
plot(density(ts(imag_parts)), pty=1, main="Posterior Im" )

for(i in 1:length(density_plots) )
{  
  idx = density_plots[i]
  mcmc_data = read.table(list_of_files[idx])
  imag_parts = mcmc_data[8:8]
  x <- head( readLines(list_of_files[idx]))
  meta_data_string <- grep( pattern=tune_variable, x, value=TRUE)
  lines(density(ts(imag_parts)), lty=i+1)
  vals[i+1] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  print(vals[i+1])
}
legend('topright',as.character(vals), bty='n', 
       cex=.75, lty=c(1, density_plots) )


#mtext("Path Length", outer=TRUE, cex=1.5)

