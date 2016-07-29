library(data.table)
par(mar=c(1,1,1,1))
par(mfrow = c(2, 2))
setwd("~/inferring_geometries/output")

lf_pathlength_var = read.table("lf_parallelpaths_variance.txt")
plot( lf_pathlength_var[c("V1","V6")] )

lf_pathlength_var = read.table("lf_pathlength_variance.txt")
plot( lf_pathlength_var[c("V2","V6")] )

lf_pathlength_var = read.table("lf_obssigma_variance.txt")
plot( lf_pathlength_var[c("V3","V6")] )

lf_pathlength_var = read.table("lf_diffsigma_variance.txt")
plot( lf_pathlength_var[c("V4","V6")] )