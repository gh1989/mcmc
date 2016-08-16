setwd("C:/cygwin64/home/Gregg/mcmc/output")

library(reshape2)
library(ggplot2)

DF = read.table("ts_smc.txt", header=FALSE)
obs = read.table("ts_smc_observation.txt", header=FALSE)

num_particles = length(names(obs))/2

plot(obs, type="p", col="blue")

num_particles = length(names(DF))/2

for(i in 1:num_particles)
  points(DF[(2*i-1):(2*i)], type="l", col=i)
