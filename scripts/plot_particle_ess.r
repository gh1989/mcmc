par(mfrow=c(1,3) )
par(mar=c(4,2,2,4))

folder = "number_particles"

size_runs = 8
the_var = "P"
the_pattern = "_path_length"
source("C:/cygwin64/home/Gregg/mcmc/scripts/ess_particle.r")

size_runs = 5
the_var = "M"
the_pattern = "_extra_data_ratio"
source("C:/cygwin64/home/Gregg/mcmc/scripts/ess_particle.r")

folder = "number_particles"

size_runs = 4
the_var = "Z"
the_pattern = "_cutoff"
source("C:/cygwin64/home/Gregg/mcmc/scripts/ess_particle_dimension.r")
