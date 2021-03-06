library(zoo)
library(mcmcse)

par(mfrow=c(5,4) ,oma = c(0, 0, 2, 0) )
par(mar=c(4,2,2,4))

folder = "LFMCMC"

first_plot = 5
density_plots = c(2,3)
ptys=c(2,3)
tune_variable = "parallel_paths"
nice_name = "K"
log_xy = 'xy'
source("C:/cygwin64/home/Gregg/mcmc/scripts/lf_mcmc_tuning_o.r")

first_plot = 4
density_plots = c(2,3)
ptys=c(2,3)
tune_variable = "path_length"
nice_name = "P"
log_xy = 'xy'
source("C:/cygwin64/home/Gregg/mcmc/scripts/lf_mcmc_tuning_o.r")

first_plot = 4
density_plots = c(5,6)
ptys=c(2,3,4,5)
tune_variable = "extra_data_ratio"
nice_name = "M"
log_xy = 'xy'
source("C:/cygwin64/home/Gregg/mcmc/scripts/lf_mcmc_tuning_o.r")

first_plot = 4
density_plots = c(3,2)
ptys=c(2,3)
tune_variable = "observation_noise_sigma"
nice_name = expression(sigma[o])
log_xy = 'xy'
source("C:/cygwin64/home/Gregg/mcmc/scripts/lf_mcmc_tuning_o.r")

first_plot = 1
density_plots = c(2,3)
ptys=c(2,3)
tune_variable = "real_sigma"
nice_name = expression(sigma[d])
log_xy = 'xy'
source("C:/cygwin64/home/Gregg/mcmc/scripts/lf_mcmc_tuning_o.r")