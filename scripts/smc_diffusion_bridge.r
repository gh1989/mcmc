# Outputs for pMCMC time series.

diffusionPlot <- function( smc_diffusion_bridge )
{
  
  plot(ts(smc_diffusion_bridge[1:1]))
  plot(ts(smc_diffusion_bridge[2:2]))
  
  acf(smc_diffusion_bridge[1:1])
  acf(smc_diffusion_bridge[2:2])
  
}
setwd("~/inferring_geometries/output")
par(mar=c(4,2,2,4))
par(mfrow = c(4, 2))
smc_observed = read.table("smc_bridge_observed.txt")
smc_diffusion_bridge = read.table("smc_bridge_timeseries.txt")
diffusionPlot(smc_diffusion_bridge)
diffusionPlot(smc_observed)

