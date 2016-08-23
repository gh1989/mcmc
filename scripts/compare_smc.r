par(mar=c(4,2,2,4))
par(mfrow = c(1, 2))

DF = read.table("ts_smc.txt", header=FALSE)
obs = read.table("ts_smc_observation.txt", header=FALSE)
probs = read.table("ts_smc_probs.txt", header=FALSE)

source("C:/cygwin64/home/Gregg/mcmc/scripts/smc_trajectories.r")

DF = read.table("ts_smc_bridge.txt", header=FALSE)
obs = read.table("ts_smc_bridge_observation.txt", header=FALSE)
probs = read.table("ts_smc_bridge_probs.txt", header=FALSE)

source("C:/cygwin64/home/Gregg/mcmc/scripts/smc_trajectories.r")

