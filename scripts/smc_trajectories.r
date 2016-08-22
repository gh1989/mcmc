
DF = read.table("ts_smc.txt", header=FALSE)
obs = read.table("ts_smc_observation.txt", header=FALSE)
probs = read.table("ts_smc_probs.txt", header=FALSE)

num_particles = length(names(obs))/2

min_x = min(obs[1:1])-0.01
max_x = max(obs[1:1])+0.075
min_y = min(obs[2:2])-0.01
max_y = max(obs[2:2])+0.025

xlim <- c( min_x, max_x )
ylim <- c( min_y, max_y )

x <- seq( min_x, max_x, len=100 )
y <- seq( min_y, max_y, len=100 )
z <- matrix( sin( 2*pi*(expand.grid(x,y)$Var1 + expand.grid(x,y)$Var2)) + cos( 2*pi*(expand.grid(x,y)$Var1 + expand.grid(x,y)$Var2)),length(x),length(y))
Lab.palette <- colorRampPalette(c("Black","Red", "Yellow"), space = "rgb")
filled.contour(x,y,z, xlab="x", ylab="y", xlim<-xlim, ylim<-ylim, nlevels=128) #, col=Lab.palette(128))

num_particles = length(names(DF))/2

max_probs = max(probs)

for(i in 1:num_particles)
{
  d = as.numeric(1 - probs[i] / max_probs )
  if ((probs[i] < 0.40*max_probs) && (probs[i]>0.001*max_probs))
  {
    p_col = rgb(0.45,0.45,0.45)
    points(DF[(2*i-1):(2*i)], type="l", col=p_col)
  }

}

for(i in 1:num_particles)
{
    d = as.numeric(1 - probs[i] / max_probs )
    if (probs[i] > 0.40*max_probs)
    {
        p_col = rgb(0,1-d,0)
        points(DF[(2*i-1):(2*i)], type="l", col=p_col, pch=3)
    }
}

points(obs[1:2], type="b", lty=3, col="blue", pch=20 )
#text(obs[1:2], as.character(seq(1,8)), col="white", cex=0.75, offset=0.25)
legend("topleft", c("Likely", "Unlikely", "Observations"), col=c("green","grey", "blue"),  lty=c(1,1,3), pch=c(3,0,20), cex=.75 )
