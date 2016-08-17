library(reshape2)
library(ggplot2)

DF = read.table("ts_smc.txt", header=FALSE)
obs = read.table("ts_smc_observation.txt", header=FALSE)
probs = read.table("ts_smc_probs.txt", header=FALSE)

num_particles = length(names(obs))/2

plot(obs[1:2], type="p", col="blue")

num_particles = length(names(DF))/2

max_probs = max(probs)

for(i in 1:num_particles)
{
  d = as.numeric(1 - probs[i] / max_probs )
  
  if (probs[i] < 0.90*max_probs)
  {
    p_col = rgb(1-d,1-d,1-d)
  }
  else
    p_col = rgb(0,1-d,0)

  if (probs[i] > 0.80*max_probs)
  {
    points(DF[(2*i-1):(2*i)], type="o", col=p_col)
  }
  
}

points(obs[1:2], type="p", col="blue")